classdef HypersphericalDiracDistribution < AbstractHypersphereSubsetDiracDistribution & AbstractHypersphericalDistribution
    % Dirac distribution on the hypersphere with dirac positions d and
    % weights w.
    
    methods   
        function p = plot(this, varargin)
            % Create an appropriate plot
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot/surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            switch this.dim
                case 2
                    % use polar coordinates, plot angle->pdf(angle)
                    p = stem(atan2(this.d(2,:), this.d(1,:)), this.w, varargin{:});
                case 3
                    % plot points on sphere
                    p = scatter3(this.d(1,:),this.d(2,:),this.d(3,:), this.w*size(this.d,2)*20, varargin{:});
                otherwise
                    error('Plotting for this dimension is currently not supported');
            end                        
        end
        
        function wd = toWD(this)
            % Convert to a WD distribution (only in 2D case)
            %
            % Returns:
            %   wd (WDDistribution)
            %       WDDistribution with same parameters
            assert(this.dim == 2);
            wd = WDDistribution(atan2(this.d(2,:), this.d(1,:)), this.w);
        end   
        
        function mu = meanDirection(this)
            vecSum = sum(this.d.*this.w, 2);
            mu = vecSum/norm(vecSum);
        end
        
        function [s,info] = sampleDeterministicLCD(this, n)
            % Spherical LCD DMD-to-DMD deterministic sampling algorithm
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            %   info (struct)
            %       more details: generated reference samples, optimization output 
            %
            % https://isas.iar.kit.edu/pdf/IFAC20_Frisch.pdf
            % Daniel Frisch, Kailai Li, and Uwe D. Hanebeck 
            % Optimal Reduction of Dirac Mixture Densities on the 2-sphere 
            % Proceedings of the 1st Virtual IFAC World Congress (IFAC-V 2020), July, 2020.  
            arguments
                this (1,1) HypersphericalDiracDistribution
                n    (1,1) double {mustBeNonnegative, mustBeReal, mustBeInteger}
            end
            assert(this.dim==3, 'Spherical LCD Distance is currently not supported in this dimension.') 
            
            % Cartesian to spherical coordinates, (3 x n) [x;y;z] --> (2 x n) [theta;phi]
            c2s = @(c) [atan2(vecnorm(c(1:2,:),2,1),c(3,:)); atan2(c(2,:),c(1,:))];
            % spherical to Cartesian coordinates, (2 x n) [theta;phi] --> (3 x n) [x;y;z]
            s2c = @(s) [sin(s(1,:)).*cos(s(2,:)); sin(s(1,:)).*sin(s(2,:)); cos(s(1,:))];
            % geodesic distance from spherical coordinates, ((2 x nx), (2 x ny)) [theta;phi] --> (nx x ny), more precise than acos(sdot(sx,xy))
            sgeo = @(sx,sy) 2 * asin(sqrt( sin((sx(1,:)'-sy(1,:))/2).^2 + sin(sx(1,:)').*sin(sy(1,:)).*sin((sx(2,:)'-sy(2,:))/2).^2 ));
            % Ei function (without the imaginary part that cancels out anyway)
            ei = @(x) real( -expint(-x) ); % TODO use mex from nonlinearestimation toolbox
            
            % obtain initial guess
            s0 = this.sample(n); % (3 x n) [x;y;z]
            % convert to spherical coordinates for optimization
            sRef_Spher = c2s(this.d); % (2 x nRef) [theta;phi]
            s0_Spher   = c2s(s0);     % (2 x n)    [theta;phi]
            % weights
            wRef = this.w;      % (1 x nRef)
            w0   = ones(1,n)/n; % (1 x n)
            
            switch('Normal') % Normal, HighPrecision, GradientCheck 
                case 'Normal'
                    options = optimoptions('fminunc', 'Display','notify-detailed', 'SpecifyObjectiveGradient',true );
                case 'HighPrecision'
                    options = optimoptions('fminunc', 'Display','notify-detailed', 'SpecifyObjectiveGradient',true, 'OptimalityTolerance',1e-12, 'MaxIterations',2000, 'StepTolerance',1e-12 );
                case 'GradientCheck'
                    options = optimoptions('fminunc', 'Display','iter', 'SpecifyObjectiveGradient',true, 'CheckGradients',true, 'FiniteDifferenceType','central', 'FiniteDifferenceStepSize',1e-4);
                otherwise
                    error('Wrong switch')
            end

            % optimization
            [sOpt_Spher,fval,exitflag,output] = fminunc(@fmindist, s0_Spher, options);
            s = s2c(sOpt_Spher);
            
            if nargout > 1
                % return additional information
                info = struct();
                info.samplesRef.Spherical = sRef_Spher;
                info.samplesRef.Cartesian = this.d;
                info.samplesInit.Spherical = s0_Spher;
                info.samplesInit.Cartesian = s0;
                info.optimization.samplesOpt.Spherical = sOpt_Spher;
                info.optimization.samplesOpt.Cartesian = s;
                info.optimization.fval = fval;
                info.optimization.exitflag = exitflag;
                info.optimization.output = output;
            end

            function [D,Diff] = fmindist(sOpt_Spher)
                % scalar distance function for fmincon or fminunc.
                %
                % Parameters:
                %   xOpt : state vector (2 x n) 
                %
                % Returns:
                %   D    : (1 x 1) distance value
                %   Diff : (2 x n) gradient ∂D/∂sOpt 
                
                % calculate distance and gradients for sample pairs
                [DXX, DXX_dX, DXX_dY] = DXYfun(sOpt_Spher,sOpt_Spher,w0,w0);   % (1 x 1), (2 x nOpt), (2 x nOpt)
                [DXY, DXY_dX,   ~   ] = DXYfun(sOpt_Spher,sRef_Spher,w0,wRef); % (1 x 1), (2 x nOpt), (2 x nRef)
                % derivative of distance w.r.t. sOpt
                DXXDiff = DXX_dX + DXX_dY; % (2 x nOpt)
                DXYDiff = DXY_dX;          % (2 x nOpt)
                % calculate overall distance and derivative 
                D = DXX - 2*DXY;  % (1 x 1)  (YY contribution not necessary as sRef is fixed)
                Diff = DXXDiff - 2*DXYDiff; % (2 x n) (size of sOpt)
            end
            
            function [DXY, DiffXY_dX, DiffXY_dY] = DXYfun(X,Y,wx,wy)
                % Weighted distance measure between individual points
                % (equally weighted samples)
                % Can be used to calculate DXX, DXY, DYY
                %
                % Parameters:
                %   X      : (2 X nX) [theta;phi]
                %   Y      : (2 X nY) [theta;phi]
                %   wx     : (1 x nX)
                %   wy     : (1 x nY)
                % Returns: 
                %   DXY    : (1 x 1 ) distance value
                %   DXY_dX : (2 x nX) [theta;phi] ∂DXY/∂X 
                %   DXY_dY : (2 X nY) [theta;phi] ∂DXY/∂Y 
                
                % calculate geodesic distance d
                geodist = sgeo(X,Y); % (nX x nY) 
                % calculate unweighted distance measure function Ds() from d
                [Ds_XY,DiffXYMat] = Ds_Mises(geodist); % (nX x nY), (nX x nY) 
                % Apply sample weights
                DXY = wx' .* wy .* Ds_XY; % (nX x nY)
                % Sum them up to total distance value
                DXY = sum(DXY(:)); % [1 x 1]
                
                % calculate analytical gradients
                % part of ∂DUV/∂XY without ∂dc/∂XY
                DiffXYMat = - wx' .* wy .* DiffXYMat ./ sin(geodist);  % (nX x nY) 
                % set parts with d=0 and 0/sin(0)=NaN to zero
                ind = geodist<=0+1e-5;
                DiffXYMat(ind) = 0;
                DiffXYMat = permute(DiffXYMat, [3,1,2]); % (1 x nX x nY)
                % ∂dc/∂XY
                [Diff_dc_X,Diff_dc_Y] = dotXY_Diff_Spher(X,Y); % (2 x nX x nY), (2 x nX x nY)  [theta;phi]
                % ∂DXY/∂X
                DiffXY_dX = sum( DiffXYMat .* Diff_dc_X, 3); % (2 x nX)      [theta;phi]
                % ∂DXY/∂Y
                DiffXY_dY = sum( DiffXYMat .* Diff_dc_Y, 2); % (2 x 1 x nY)  [theta;phi]
                DiffXY_dY = DiffXY_dY(:,:);   % (2 x nY)  [theta;phi] 
            end
            
            function [DiffX,DiffY] = dotXY_Diff_Spher(X,Y)
                % Analytical gradient of dot product between two sets of points on sphere in Spherical coordinates
                %
                % Parameters:
                %   X : (2 x nX) [theta;phi]
                %   Y : (2 x nY) [theta;phi]
                %
                % Returns:
                %   [DiffX,DiffY] : (2 x nX x nY), (2 x nX x nY) [theta;phi]
                
                assert(size(X,1)==2)
                % variables theta and phi
                x_th  = X(1,:);                   % (1 x nX x  1)
                y_th  = permute(Y(1,:), [1,3,2]); % (1 x  1 x nY) 
                x_phi = X(2,:);                   % (1 x nX x  1)
                y_phi = permute(Y(2,:), [1,3,2]); % (1 x  1 x nY)
                % necessary sine and cosine terms
                cos_x_th = cos(x_th);
                sin_x_th = sin(x_th);
                cos_y_th = cos(y_th);
                sin_y_th = sin(y_th);
                cos_xyphi = cos(x_phi-y_phi);
                sin_xyphi = sin(x_phi-y_phi);
                % definition of D for reference
                %      D  =   sin(x_th) .* sin(y_th) .* cos(x_phi-y_phi) + cos(x_th).*cos(y_th);
                % diff of D w.r.t. x_th and y_th
                DiffX_th  =   cos_x_th .* sin_y_th .* cos_xyphi - sin_x_th.*cos_y_th; % (1 x nX x nY)
                DiffY_th  =   sin_x_th .* cos_y_th .* cos_xyphi - cos_x_th.*sin_y_th; % (1 x nX x nY)
                % diff w.r.t. x_phi and y_phi
                DiffX_phi = - sin_x_th .* sin_y_th .* sin_xyphi; % (1 x nX x nY)
                DiffY_phi =   sin_x_th .* sin_y_th .* sin_xyphi; % (1 x nX x nY)
                % combine differences w.r.t. th and phi to array
                DiffX = cat(1, DiffX_th, DiffX_phi); % (2 x nX x nY)
                DiffY = cat(1, DiffY_th, DiffY_phi); % (2 x nX x nY)
            end
            
            function [D,Diff] = Ds_Mises(d)
                % Ds(d) function for Von-Mises Fisher Kernel
                % Calculates Distance between two samples of geodesic distance d 
                %
                % Parameters:
                %   d      : (n1 x n2) 
                % Returns: 
                %   D      : (n1 x n2) D 
                %   Diff   : (n1 x n2) ∂D/∂d
                 
                % integral_(b1=0.001)^(b2=Inf)
                % b2 -> Inf, then the indefinite integral real part turns to zero: 
                D2 = 0;  Diff2 = 0;
                b1 = .001; 
                [D1,Diff1] = Ds_Mises_indefinite(d,b1);
                D = D2-D1;
                Diff = Diff2-Diff1;
                
                function [D,Diff] = Ds_Mises_indefinite(d,b)
                    % define intermediate variables
                    cosd2 = cos(d/2);
                    c1 = +2 * (cosd2 - 1);
                    c2 = -2 * (cosd2 + 1);
                    Ei_c1 = ei(c1.*b);
                    Ei_c2 = ei(c2.*b);
                    % case for d!=0:
                    D = pi./cosd2 .* ( c1.*Ei_c1 - c2.*Ei_c2  +  1./b.*(exp(c2.*b)-exp(c1.*b)) );
                    Diff = tan(d/2)/2 .* D  -  pi .* tan(d/2) .* ( Ei_c1 + Ei_c2 );
                    % case for d=0 (equivalent to c1=0):
                    ind = d<=0+1e-5;
                    D(ind) = pi * ( 4*ei(-4*b) + 1./b.*(exp(-4*b)-1) );
                    Diff(ind) = 0;
                    % case for d=pi (equivalent to c2=0):
                    ind = d>=pi-1e-5;
                    D(ind) = 4*pi*ei(-2*b);
                    Diff(ind) = 0;
                end
            end
        end
    end
    methods (Static)
        function f = fromDistribution(distribution, noOfSamples)
            arguments
                distribution (1,1) AbstractHypersphericalDistribution
                noOfSamples (1,1) {mustBePositive, mustBeInteger}
            end
            f = HypersphericalDiracDistribution(...
                distribution.sample(noOfSamples),ones(1,noOfSamples)/noOfSamples);
        end
    end
end

