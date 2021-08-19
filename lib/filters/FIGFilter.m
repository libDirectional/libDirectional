classdef FIGFilter < AbstractCircularFilter & AbstractGridFilter
    % Implementation of the Fourier-interpreted grid filter.
    % see Florian Pfaff, Kailai Li, and Uwe D. Hanebeck,
    % Fourier Filters, Grid Filters, and the Fourier-Interpreted Grid Filter,
    % Proceedings of the 22nd International Conference on Information Fusion (Fusion 2019), Ottawa, Canada, July, 2019.
    
    methods
        function this = FIGFilter(noOfCoefficients, enforcePdfNonnegative)
            arguments
                noOfCoefficients {mustBeInteger,mustBePositive}
                enforcePdfNonnegative logical = true
            end
            this.gd = FIGDistribution.fromDistribution(CircularUniformDistribution(), noOfCoefficients, enforcePdfNonnegative);
        end
        
        function predictIdentity(this, dSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by dSys.
            %
            % Parameters:
            %   dSys (AbstractCircularDistribution or samples on grid)
            %       distribution of additive noise or samples on a grid. If
            %       values on a grid are given, take care that they are
            %       not rooted yet. Take care that if
            %       linspace(0,2*pi,noCoeffs) is used, that the last grid
            %       point has to be discarded. For example
            %       xvals=0:2*pi/noCoeffs:2*pi-2*pi/noCoeffs;
            assert(isa(dSys, 'AbstractCircularDistribution'))
            % If a distribution on a circle is given, use Hadamard
            % product approach
            if ~isa(dSys, 'FIGDistribution')
                warning('Predict:automaticConversion', ...
                    'dSys is not a FIGDistribution. Transforming with a number of coefficients that is equal to that of the filter. For non-varying noises, transforming once is much more efficient and should be preferred.');
                dSys = FIGDistribution.fromDistribution(dSys, numel(this.gd.gridValues), this.gd.enforcePdfNonnegative);
            end
            this.gd = this.gd.convolve(dSys);
        end
        
        function updateIdentity(this, dMeas, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by dMeas.
            %
            % Parameters:
            %   dMeas (AbstractCircularDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            arguments
                this (1,1) FIGFilter
                dMeas (1,1) AbstractCircularDistribution
                z (1,1) double = 0
            end
            if ~isa(dMeas, 'FIGDistribution') % Include shift in transformation
                dMeas = FIGDistribution.fromFunction(@(x)dMeas.pdf(x-z), numel(this.gd.gridValues), this.gd.enforcePdfNonnegative);
            elseif z > 1e-8 % .shift is rather expensive, avoid if possible
                dMeas = dMeas.shift(z);
            end
            this.gd = this.gd.multiply(dMeas);
        end
        
        function predictNonlinear(this, f, noiseDistribution, truncateJointSqrt)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by noiseDistribution.
            % Using predictNonlinearViaTransitionDensity with a
            % pretransformed fTrans is to be preferred as the function is
            % always performed anew when using predictNonlinear.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi).
            %       fuction MUST support row vectors of different angles
            %       and output a row vector!
            %   noiseDistribution (AbstractCircularDistribution)
            %       distribution of additive noise
            
            assert(isa(noiseDistribution, 'AbstractCircularDistribution'));
            assert(isa(f, 'function_handle'));
            if nargin == 3, truncateJointSqrt = true; end
            fTrans = @(xkk, xk)reshape(noiseDistribution.pdf(xkk(:)'-f(xk(:)')), size(xk));
            this.predictNonlinearViaTransitionDensity(fTrans, truncateJointSqrt);
            
        end
        
        function predictNonlinearViaTransitionDensity(this, fTrans, truncateJointSqrt) %#ok<INUSD>
            error('Currently unsupported');
        end
    end
    
end
