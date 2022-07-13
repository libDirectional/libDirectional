classdef (Abstract) AbstractHypersphericalDistribution < AbstractHypersphereSubsetDistribution
    % Abstract base class for distributions on the hypershere (S^dim)
    % Convention for Dimension: dim=2 is a circle, dim=3 a sphere.
    
    methods
        function h = plot(this, faces, gridFaces)
            % Plots the pdf of a hyperspherical distribution
            %
            % Parameters:
            %   faces (scalar)
            %       Number of faces for 3D Plot (default 100).
            %   gridFaces (scalar)
            %       Number of grid faces for 3D Plot (default 20, 0 to disable).
            arguments
                this (1,1) AbstractHypersphericalDistribution
                faces (1,1) double {mustBeNonnegative,mustBeInteger} = 100
                gridFaces (1,1) double {mustBeNonnegative,mustBeInteger} = 20
            end
            switch this.dim
                case 2
                    % use polar coordinates, plot angle->pdf(angle)
                    phi = linspace(0,2*pi,320);
                    x = [cos(phi); sin(phi)];
                    p = this.pdf(x);
                    h = plot(phi, p);
                case 3
                    % plot sphere, pdf is represented by color on sphere
                    % generate spheres
                    [xSphereOuter, ySphereOuter, zSphereOuter]= sphere(gridFaces);
                    [xSphereInner, ySphereInner, zSphereInner]= sphere(faces);
                    
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
                this (1,1) AbstractHypersphericalDistribution
            end
            if this.dim<=4
                mu = meanDirectionNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi;pi*ones(this.dim-2,1)]]);
            else
                % use monte carlo integration
                Sd = this.getManifoldSize();
                
                n = 10000; % number of samples for integration
                r = HypersphericalUniformDistribution(this.dim).sample(n);
                p = this.pdf(r);
                
                mu = r*p'/n * Sd;
            end
            if norm(mu)<1e-9
                warning('Density may not have actually have a mean direction because integral yields a point very close to the origin.')
            end
            mu = mu/norm(mu);
        end

        function m = momentNumerical(this)
            m = momentNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi;pi*ones(this.dim-2,1)]]);
        end
        
        function i = integralNumerical(this)
            if this.dim<=4
                i = integralNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi;pi*ones(this.dim-2,1)]]);
            else
                % use monte carlo integration
                n = 10000; % number of samples for integration
                r = HypersphericalUniformDistribution(this.dim).sample(n);
                p = this.pdf(r);
                Sd = AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim);
                i = sum(p)/n * Sd; % average value of pdf times surface area of unit sphere
            end
        end
        
        function i = entropyNumerical(this)
            arguments
                this (1,1) AbstractHypersphericalDistribution
            end
            i = entropyNumerical@AbstractHypersphereSubsetDistribution(this, [zeros(this.dim-1,1), [2*pi;pi*ones(this.dim-2,1)]]);
        end
        
        function m = modeNumerical(this)
            % Calculate the mode by finding maximum of PDF numerically 
            % Returns:
            %   m (column vector)
            %       mode of the distribution
            arguments
                this (1,1) AbstractHypersphericalDistribution
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
                this (1,1) AbstractHypersphericalDistribution
                other (1,1) AbstractHypersphericalDistribution
            end
            dist = hellingerDistanceNumerical@AbstractHypersphereSubsetDistribution(this, other, [zeros(this.dim-1,1), [2*pi;pi*ones(this.dim-2,1)]]);
        end

        function dist = totalVariationDistanceNumerical(this, other)
            arguments
                this (1,1) AbstractHypersphericalDistribution
                other (1,1) AbstractHypersphericalDistribution
            end
            dist = totalVariationDistanceNumerical@AbstractHypersphereSubsetDistribution(this, other, [zeros(this.dim-1,1), [2*pi;pi*ones(this.dim-2,1)]]);
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
                startPoint (:,1) double = this.meanDirection()
                burnIn (1,1) double = 10
                skipping (1,1) double = 5
            end
            if nargin(proposal)==0
                normalize = @(x) x/norm(x);
                % Think about this carefully. In general, not good for
                % multimodality.
                %proposal = @(x) normalize(x + mvnrnd(zeros(this.dim,1),eye(this.dim))'); 
                proposal = @(x) normalize(x + randn(this.dim,1)); 
            end
            s = sampleMetropolisHastings@AbstractDistribution(this, n, proposal, startPoint, burnIn, skipping);
        end
        
        function s = sampleDeterministic(this)
            % Deterministic sampling
            % Fall back to spherical LCD DMD-to-DMD deterministic sampling by default 
            % Use small number of samples according to UKF
            %
            % Returns:
            %   s (dim x n)
            %       one sample per column
            n = (this.dim-1)*2 + 1;
            s = sampleDeterministicLCD(this, n);
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
            % https://www.overleaf.com/read/bjchtykkytdf
            % Daniel Frisch, Kailai Li, and Uwe D. Hanebeck 
            % Optimal Reduction of Dirac Mixture Densities on the 2-sphere 
            % IFAC 2020, Berlin 
            arguments
                this (1,1) AbstractHypersphericalDistribution
                n    (1,1) double {mustBeNonnegative, mustBeReal, mustBeInteger}
            end
            assert(this.dim==3, 'Spherical LCD Distance is currently not supported in this dimension.') 
            % obtain stochastic reference samples
            nRef = n*100; 
            sRef = this.sample(nRef); % (3 x nRef) [x;y;z]
            % DMD object, with equally weighted samples 
            DMD = HypersphericalDiracDistribution(sRef);
            % get LCD sample reduction
            [s,info] = DMD.sampleDeterministicLCD(n);
        end
        
        function s = getManifoldSize(this)
            arguments
                this (1,1) AbstractHypersphericalDistribution
            end
            s = AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim);
        end
    end
    
    methods (Static)
        function surfaceArea = computeUnitSphereSurface(dimension)
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
            
            % Use switch to avoid calling gamma function for lower
            % dimensions
            switch dimension
                case 2
                    surfaceArea = 2 * pi;
                case 3
                    surfaceArea = 4 * pi;
                case 4
                    surfaceArea = 2 * pi^2;
                otherwise
                    surfaceArea = 2 * pi^((dimension)/2)/gamma((dimension)/2); 
            end
        end
        
        function h = plotSphere(varargin)
            % Plots a sphere
            %
            % Returns:
            %   h (handle)
            %       plot handle
            [x,y,z]=sphere(150); % Create smooth sphere
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

