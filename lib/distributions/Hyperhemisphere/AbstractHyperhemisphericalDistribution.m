classdef (Abstract) AbstractHyperhemisphericalDistribution < AbstractHypersphereSubsetDistribution
    % Abstract base class for distributions on the hypershere (S^dim)
    
    methods
        function m = mean(this)
            m = this.meanAxis();
        end
        function h = plot(this, faces, gridFaces)
            % Plots the pdf of a hyperspherical distribution
            %
            % Parameters:
            %   faces (scalar)
            %       Number of faces for 3D Plot (default 100).
            %   gridFaces (scalar)
            %       Number of grid faces for 3D Plot (default 20, 0 to disable).
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
                faces (1,1) double {mustBeNonnegative,mustBeInteger} = 100
                gridFaces (1,1) double {mustBeNonnegative,mustBeInteger} = 20
            end
            switch this.dim
                case 2
                    % use polar coordinates, plot angle->pdf(angle)
                    phi = linspace(0,pi,320);
                    x = [cos(phi); sin(phi)];
                    p = this.pdf(x);
                    h = plot(phi, p);
                case 3
                    % plot sphere, pdf is represented by color on sphere
                    % generate spheres
                    [xSphereOuter, ySphereOuter, zSphereOuter] = sphere(gridFaces);
                    xSphereOuter=reshape(xSphereOuter(zSphereOuter(:)>=0),[gridFaces/2+1,gridFaces+1]);
                    ySphereOuter=reshape(ySphereOuter(zSphereOuter(:)>=0),[gridFaces/2+1,gridFaces+1]);
                    zSphereOuter=reshape(zSphereOuter(zSphereOuter(:)>=0),[gridFaces/2+1,gridFaces+1]);
                    
                    [xSphereInner, ySphereInner, zSphereInner] = sphere(faces);
                    xSphereInner=reshape(xSphereInner(zSphereInner(:)>=0),[faces/2+1,faces+1]);
                    ySphereInner=reshape(ySphereInner(zSphereInner(:)>=0),[faces/2+1,faces+1]);
                    zSphereInner=reshape(zSphereInner(zSphereInner(:)>=0),[faces/2+1,faces+1]);
                    
                    
                    % evaluate p.d.f.
                    cSphere=reshape(this.pdf([xSphereInner(:)';ySphereInner(:)';zSphereInner(:)']),size(xSphereInner));
                    
                    % resize inner sphere
                    xSphereInner = 0.99 * xSphereInner;
                    ySphereInner = 0.99 * ySphereInner;
                    zSphereInner = 0.99 * zSphereInner;
                    
                    holdStatus=ishold;
                    % Plot spheres
                    h = [];
                    if gridFaces > 0
                        h = [h, ...
                            surf(xSphereOuter, ySphereOuter, zSphereOuter, max(cSphere,[],[1,2])*ones(size(xSphereOuter)), 'FaceColor', 'none')];
                        hold on;
                    end
                    h = [h, surf(xSphereInner, ySphereInner, zSphereInner, cSphere,'EdgeColor', 'none')];
                    axis equal
                    colorbar
                    if ~holdStatus
                        hold off
                    end
                otherwise
                    error('Cannot plot hyperspherical distribution with this number of dimensions');
            end
        end
        
        function mu = meanDirectionNumerical(this)
            % Calculate mean direction of pdf
            % Returns:
            %   mu (vector)
            %       mean direction
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
            end
            warning(['The result is the mean direction on the upper hemisphere ',...
            'along the last dimension. It is not a mean of a symmetric ',...
            'distribution, which would not have a proper mean. ',...
            'It is also not one of the modes of the symmetric distribution ',...
            'since it is biased toward [0;...;0;1] because the lower half ',...
            'is considered inexistant.']);
            if this.dim==2
                mu = meanDirectionNumerical@AbstractHypersphereSubsetDistribution(this, [0, pi]);
            elseif this.dim<=4
                mu = meanDirectionNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi; pi*ones(this.dim-3,1);pi/2]]); % 2*pi, pi, ..., pi/2
            else
                % use monte carlo integration
                Sd = this.getManifoldSize();
                
                n = 10000; % number of samples for integration
                r = HyperhemisphericalUniformDistribution(this.dim).sample(n);
                p = this.pdf(r);
                
                mu = r*p'/n * Sd;
            end
            if norm(mu)<1e-9
                warning('Density may not have actually have a mean direction because integral yields a point very close to the origin.')
            end
            mu = mu/norm(mu);
        end

        function m = momentNumerical(this)
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
            end
            if this.dim==2
                m = momentNumerical@AbstractHypersphereSubsetDistribution(this, [0, pi]);
            else
                m = momentNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi; pi*ones(this.dim-3,1);pi/2]]);
            end
        end
        
        function i = integralNumerical(this)
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
            end
            if this.dim==2
                i = integralNumerical@AbstractHypersphereSubsetDistribution(this, [0, pi]);
            elseif this.dim<=4
                i = integralNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi; pi*ones(this.dim-3,1);pi/2]]);
            else
                % use monte carlo integration
                n = 10000; % number of samples for integration
                r = HyperhemisphericalUniformDistribution(this.dim).sample(n);
                p = this.pdf(r);
                Sd = AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim);
                i = sum(p)/n * Sd; % average value of pdf times surface area of unit sphere
            end
        end
        
        function i = entropyNumerical(this)
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
            end
            if this.dim==2
                i = entropyNumerical@AbstractHypersphereSubsetDistribution(this, [0, pi]);
            else
                i = entropyNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi; pi*ones(this.dim-3,1);pi/2]]);
            end
        end
        
        function m = modeNumerical(this)
            % Calculate the mode by finding maximum of PDF numerically 
            % Returns:
            %   m (column vector)
            %       mode of the distribution 
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
            end
            fun = @(s) -this.pdf(polar2cart(s)); % objective function 
            s0 = rand(this.dim-1,1)*pi; % initial point 
            options = optimoptions('fminunc', 'Display','notify-detailed', 'OptimalityTolerance',1e-12, 'MaxIterations',2000, 'StepTolerance',1e-12 ); 
            % find maximum of density by numerical optimization in spherical coordinates 
            m = fminunc(fun, s0, options);  
            m = polar2cart(m); % convert optimum to Cartesian coordinates 
        end

        function dist = hellingerDistanceNumerical(this, other)
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
                other (1,1) AbstractHyperhemisphericalDistribution
            end
            if this.dim==2
                dist = hellingerDistanceNumerical@AbstractHypersphereSubsetDistribution(this, other, [0, pi]);
            else
                dist = hellingerDistanceNumerical@AbstractHypersphereSubsetDistribution(this, other, [zeros(this.dim-1,1), [2*pi; pi*ones(this.dim-3,1);pi/2]]);
            end
        end

        function dist = totalVariationDistanceNumerical(this, other)
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
                other (1,1) AbstractHyperhemisphericalDistribution
            end
            if this.dim==2
                dist = totalVariationDistanceNumerical@AbstractHypersphereSubsetDistribution(this, other, [0, pi]);
            else
                dist = totalVariationDistanceNumerical@AbstractHypersphereSubsetDistribution(this, other, [zeros(this.dim-1,1), [2*pi; pi*ones(this.dim-3,1);pi/2]]);
            end
        end
        
        function s = sampleMetropolisHastings(this, n, proposal, startPoint, burnIn, skipping)
            % Metropolis Hastings sampling algorithm
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            %
            % Hastings, W. K. 
            % Monte Carlo Sampling Methods Using Markov Chains and Their Applications 
            % Biometrika, 1970, 57, 97-109
            arguments
                this (1,1) AbstractDistribution
                n (1,1) {mustBePositive,mustBeInteger}
                proposal (1,1) function_handle = @()[] % Will be replaced by default below
                startPoint (:,1) double = this.mode()
                burnIn (1,1) double = 10
                skipping (1,1) double = 5
            end
            if nargin(proposal)==0
                normalize = @(x) x/norm(x);
                % This is just an approximation and we actually allow
                % points and their antipode to be included when the last
                % entry is zero on the last dimension. We ignore this
                % because it is generally unlikely, but changes may be
                % required if that that cause problems.
                toUpperHemisphere = @(s) (1-2*(s(end,:)<0)).*s;
                proposal = @(x) toUpperHemisphere(normalize(x + normrnd(0,1,this.dim,1))); 
            end
            s = sampleMetropolisHastings@AbstractDistribution(this, n, proposal, startPoint, burnIn, skipping);
        end
        
        function s = getManifoldSize(this)
            arguments
                this (1,1) AbstractHyperhemisphericalDistribution
            end
            s = AbstractHyperhemisphericalDistribution.computeUnitHyperhemisphereSurface(this.dim);
        end
    end
    
    methods (Static)
        function surfaceArea = computeUnitHyperhemisphereSurface(dimension)
            % Computes surface area of (d-1)-sphere
            % see http://en.wikipedia.org/wiki/N-sphere#Volume_and_surface_area
            % dimension = 2 => circle
            % dimension = 3 => sphere
            % ...
            % 
            % Parameters:
            %   dimension (scalar)
            %       the dimension of the sphere
            % Returns:
            %   surfaceArea (scalar)
            %       the surface area of the sphere
            surfaceArea = 0.5 * AbstractHypersphericalDistribution.computeUnitSphereSurface(dimension);
        end
        
        function h = plotHemisphere(varargin)
            % Plots a sphere
            %
            % Returns:
            %   h (handle)
            %       plot handle
            resolution = 150;
            [x,y,z]=sphere(resolution); % Create smooth sphere
            
            x=reshape(x(z(:)>=0),[resolution/2+1,resolution+1]);
            y=reshape(y(z(:)>=0),[resolution/2+1,resolution+1]);
            z=reshape(z(z(:)>=0),[resolution/2+1,resolution+1]);
            
            h=mesh(x,y,z);
            skipx=10;
            skipy=10;
            x=get(h,'xdata'); % Get lines from smooth sphere
            y=get(h,'ydata');
            z=get(h,'zdata');
            delete(h)
            xKeep=x(1:skipx:end,:); % Only plot some of the lines as grid would be too fine otherwise
            yKeep=y(1:skipx:end,:);
            zKeep=z(1:skipx:end,:);
            lineHandles=line(xKeep',yKeep',zKeep','color',0.7*[1 1 1]);

            xKeep=x(:,1:skipy:end);
            yKeep=y(:,1:skipy:end);
            zKeep=z(:,1:skipy:end);
            lineHandles=[lineHandles;...
                line(xKeep,yKeep,zKeep,'color',0.7*[1 1 1])];
            for hCurr=lineHandles'
                hasbehavior(hCurr,'legend',false); % Prevent legend entry for all lines
            end

            axis equal
        end
    end 
end

