classdef (Abstract) AbstractHypersphericalDistribution < AbstractDistribution
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
            switch this.dim
                case 2
                    % use polar coordinates, plot angle->pdf(angle)
                    phi = linspace(0,2*pi,320);
                    x = [cos(phi); sin(phi)];
                    p = this.pdf(x);
                    h = plot(phi, p);
                case 3
                    % plot sphere, pdf is represented by color on sphere
                    
                    if nargin < 2
                        faces = 100;
                    end
                    
                    if nargin < 3
                        gridFaces = 20;
                    end
                    
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
                            surf(xSphereOuter, ySphereOuter, zSphereOuter, max(max(cSphere))*ones(size(xSphereOuter)), 'FaceColor', 'none')];
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
        
        function mu = meanDirection(this)
            % Calculate mean direction of pdf
            % Returns:
            %   mu (vector)
            %       mean direction
            mu = meanDirectionNumerical(this);
        end
        
        function mu = meanDirectionNumerical(this)
            % Calculate mean direction of pdf
            % Returns:
            %   mu (vector)
            %       mean direction
            mu = NaN(this.dim,1);
            switch this.dim
                case 2
                    for i=1:2
                        f = @(x) x(i,:).*this.pdf(x);
                        fAngles = @(phi) reshape(f([cos(phi(:)'); sin(phi(:)')]),size(phi));
                        mu(i) = integral(fAngles,0,2*pi, 'AbsTol', 0.01);
                    end
                case 3
                    for i=1:3
                        f = @(x) x(i,:).*this.pdf(x);
                        r = 1;

                        % spherical coordinates
                        fangles = @(phi1,phi2) f([ ...
                        r.*sin(phi1).*sin(phi2); ...
                        r.*cos(phi1).*sin(phi2); ...
                        r.*cos(phi2); ...
                        ]);

                        g = @(phi1,phi2) reshape(fangles(phi1(:)',phi2(:)').*sin(phi2(:)'),size(phi1)); % volume correcting term
                        mu(i) = integral2(g, 0, 2*pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                    end
                case 4
                    % use matlab integration
                    for i=1:4
                        f = @(x) x(i,:).*this.pdf(x);
                        r = 1;

                        % hyperspherical coordinates
                        fangles = @(phi1,phi2,phi3) f([ ...
                        r.*sin(phi1).*sin(phi2).*sin(phi3); ...
                        r.*cos(phi1).*sin(phi2).*sin(phi3); ...
                        r.*cos(phi2).*sin(phi3); ...
                        r.*cos(phi3)
                        ]);

                        g = @(phi1,phi2,phi3) fangles(phi1,phi2,phi3) .* sin(phi2).*(sin(phi3)).^2; % volume correcting term
                        ga = @(phi1,phi2,phi3) reshape(g(phi1(:)', phi2(:)', phi3(:)'), size(phi1));

                        mu(i) = integral3(ga, 0, 2*pi, 0, pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                    end
                otherwise
                    % use monte carlo integration
                    Sd = AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim);
                    
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
        
        function i = integral(this)
            % Calculate integral of pdf to check normalization
            % (should always be 1)
            % Returns:
            %   i (scalar)
            %       integral over hypersphere surface of pdf (uses
            %       approximation, not very accurate for higher dimensions)
            i = this.integralNumerical();
        end
        
        function i = integralNumerical(this)
            % Calculate integral of pdf to check normalization
            % (should always be 1)
            % Returns:
            %   i (scalar)
            %       integral over hypersphere surface of pdf (uses
            %       approximation, not very accurate for higher dimensions)
            if this.dim==2
                %use matlab integration
                f = @(phi) this.pdf([cos(phi); sin(phi)]);
                i = integral(f,0,2*pi, 'AbsTol', 0.01);
            elseif this.dim ==3
                % use matlab integration
                f = @(x) this.pdf(x);
                r=1;
                
                % spherical coordinates
                fangles = @(phi1,phi2) f([ ...
                r.*sin(phi1).*sin(phi2); ...
                r.*cos(phi1).*sin(phi2); ...
                r.*cos(phi2); ...
                ]);

                g = @(phi1,phi2) reshape(fangles(phi1(:)',phi2(:)').*sin(phi2(:)'),size(phi1)); % volume correcting term

                i = integral2(g, 0, 2*pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
            elseif this.dim==4     
                % use matlab integration
                f = @(x) this.pdf(x);
                r=1;
                
                % hyperspherical coordinates
                fangles = @(phi1,phi2,phi3) f([ ...
                r.*sin(phi1).*sin(phi2).*sin(phi3); ...
                r.*cos(phi1).*sin(phi2).*sin(phi3); ...
                r.*cos(phi2).*sin(phi3); ...
                r.*cos(phi3)
                ]);

                g = @(phi1,phi2,phi3) fangles(phi1,phi2,phi3) .* sin(phi2).*(sin(phi3)).^2; % volume correcting term
                ga = @(phi1,phi2,phi3) reshape(g(phi1(:)', phi2(:)', phi3(:)'), size(phi1));

                i = integral3(ga, 0, 2*pi, 0, pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
            else
                % use monte carlo integration
                n = 10000; % number of samples for integration
                r = HypersphericalUniformDistribution(this.dim).sample(n);
                p = this.pdf(r);
                Sd = AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim);
                i = sum(p)/n * Sd; % average value of pdf times surface area of unit sphere
            end
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically if possible, 
            % fall back to numerical calculation by default
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            result = this.entropyNumerical();
        end        
        
        function i = entropyNumerical(this)
            % Calculates the entropy numerically
            %
            % Returns:
            %   i (scalar)
            %       entropy of the distribution
            if this.dim==2
                %use matlab integration
                f = @(phi) this.pdf([cos(phi); sin(phi)]).*log(this.pdf([cos(phi); sin(phi)]));
                i = -integral(f,0,2*pi, 'AbsTol', 0.01);
            elseif this.dim ==3
                % use matlab integration
                f = @(x) this.pdf(x).*log(this.pdf(x));
                r=1;
                
                % spherical coordinates
                fangles = @(phi1,phi2) f([ ...
                r*sin(phi1).*sin(phi2); ...
                r*cos(phi1).*sin(phi2); ...
                r*cos(phi2); ...
                ]);

                g = @(phi1,phi2) fangles(phi1,phi2) .* sin(phi2); % volume correcting term
                ga = @(phi1,phi2) reshape(g(phi1(:)', phi2(:)'), size(phi1));

                i = -integral2(ga, 0, 2*pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
            elseif this.dim==4     
                % use matlab integration
                f = @(x) this.pdf(x).*log(this.pdf(x));
                r=1;
                
                % hyperspherical coordinates
                fangles = @(phi1,phi2,phi3) f([ ...
                r*sin(phi1).*sin(phi2).*sin(phi3); ...
                r*cos(phi1).*sin(phi2).*sin(phi3); ...
                r*cos(phi2).*sin(phi3); ...
                r*cos(phi3)
                ]);

                g = @(phi1,phi2,phi3) fangles(phi1,phi2,phi3) .* sin(phi2).*(sin(phi3)).^2; % volume correcting term
                ga = @(phi1,phi2,phi3) reshape(g(phi1(:)', phi2(:)', phi3(:)'), size(phi1));

                i = -integral3(ga, 0, 2*pi, 0, pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
            else
                error('not supported')
            end
        end
        
        function m = mode(this)
            % Calculate the mode by finding maximum of PDF numerically 
            % Returns:
            %   m (column vector)
            %       mode of the distribution 
            m = this.modeNumerical(); % fallback to numerical calculation (nonlinear optimization)
        end
        
        function m = modeNumerical(this)
            % Calculate the mode by finding maximum of PDF numerically 
            % Returns:
            %   m (column vector)
            %       mode of the distribution 
            fun = @(s) -this.pdf(polar2cart(s)); % objective function 
            s0 = rand(this.dim-1,1)*pi; % initial point 
            options = optimoptions('fminunc', 'Display','notify-detailed', 'OptimalityTolerance',1e-12, 'MaxIterations',2000, 'StepTolerance',1e-12 ); 
            % find maximum of density by numerical optimization in spherical coordinates 
            m = fminunc(fun, s0, options);  
            m = polar2cart(m); % convert optimum to Cartesian coordinates 
        end
        
        function s = sample(this, n)
            % Stochastic sampling
            % Fall back to Metropolis Hastings by default
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column            
            s = sampleMetropolisHastings(this, n);
        end
        
        function s = sampleMetropolisHastings(this, n)
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
            assert(isscalar(n));
            assert(n>0);
            
            burnin = 10;
            skipping = 5;
            
            totalSamples = burnin+n*skipping;
            s = zeros(this.dim,totalSamples);
            x = this.mode; 
            % A better proposal distribution could be obtained by roughly estimating
            % the uncertainty of the true distribution.
            normalize = @(x) x/norm(x);
            %proposal = @(x) normalize(x + mvnrnd(zeros(this.dim,1),eye(this.dim))'); 
            proposal = @(x) normalize(x + normrnd(0,1,this.dim,1)); 
            i=1;
            pdfx = this.pdf(x);
            while i<=totalSamples
                xNew = proposal(x); %generate new sample
                pdfxNew = this.pdf(xNew);
                a = pdfxNew/pdfx;
                if a>1
                    %keep sample
                    s(:,i)=xNew;
                    x = xNew;
                    pdfx = pdfxNew;
                    i=i+1;
                else
                    r = rand(1);
                    if a > r
                        %keep sample
                        s(:,i)=xNew;
                        x = xNew;
                        pdfx = pdfxNew;
                        i=i+1;
                    else
                        %reject sample
                    end
                end
            end
            %todo handle bimodality
            s = s(:,burnin+1:skipping:end);
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
        
        function dist = hellingerDistanceNumerical(this, other)
            % Numerically calculates the Hellinger distance to another
            % distribution.
            %
            % Parameters:
            %   other (AbstractHypersphericalDistribution)
            %       distribution to compare to
            % Returns:
            %   dist (scalar)
            %       hellinger distance of this distribution to other distribution
            assert(isa(other, 'AbstractHypersphericalDistribution'));
            assert(this.dim==other.dim,'Cannot compare distributions with different number of dimensions');
            
            % Implementation is always performed using matlab integration
            % (Implementation is similar to .integral)
            switch this.dim
                case 2
                    f = @(phi) (sqrt(this.pdf([cos(phi); sin(phi)]))-sqrt(other.pdf([cos(phi); sin(phi)]))).^2;
                    dist = 0.5*integral(f,0,2*pi, 'AbsTol', 0.01);
                case 3
                    f = @(x) (sqrt(this.pdf(x))-sqrt(other.pdf(x))).^2;
                    r=1;

                    % spherical coordinates
                    fangles = @(phi1,phi2) f([ ...
                    r.*sin(phi1).*sin(phi2); ...
                    r.*cos(phi1).*sin(phi2); ...
                    r.*cos(phi2); ...
                    ]);

                    g = @(phi1,phi2) reshape(fangles(phi1(:)',phi2(:)').*sin(phi2(:)'),size(phi1)); % volume correcting term

                    dist = 0.5*integral2(g, 0, 2*pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                case 4
                    f = @(x) (sqrt(this.pdf(x))-sqrt(other.pdf(x))).^2;
                    r=1;

                    % hyperspherical coordinates
                    fangles = @(phi1,phi2,phi3) f([ ...
                    r.*sin(phi1).*sin(phi2).*sin(phi3); ...
                    r.*cos(phi1).*sin(phi2).*sin(phi3); ...
                    r.*cos(phi2).*sin(phi3); ...
                    r.*cos(phi3)
                    ]);

                    g = @(phi1,phi2,phi3) fangles(phi1,phi2,phi3) .* sin(phi2).*(sin(phi3)).^2; % volume correcting term
                    ga = @(phi1,phi2,phi3) reshape(g(phi1(:)', phi2(:)', phi3(:)'), size(phi1));

                    dist = 0.5*integral3(ga, 0, 2*pi, 0, pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                otherwise
                    error('Numerical calculation of Hellinger distance is currently not supported for this dimension.')
            end
        end
        
        function dist = totalVariationDistanceNumerical(this, other)
            % Numerically calculates the total varation distance to another distribution
            %
            % Parameters:
            %   other (AbstractHypersphericalDistribution)
            %       distribution to compare with
            % Returns:
            %   dist (scalar)
            %       total variation distance of this distribution to other distribution
            assert(isa(other, 'AbstractHypersphericalDistribution'));
            assert(this.dim==other.dim, 'Cannot compare distributions with different number of dimensions');
            
            switch this.dim
                case 2
                    f = @(phi) abs(this.pdf([cos(phi); sin(phi)])-other.pdf([cos(phi); sin(phi)]));
                    dist = 0.5*integral(f,0,2*pi, 'AbsTol', 0.01);
                case 3
                    f = @(x) abs(this.pdf(x)-other.pdf(x));
                    r=1;

                    % spherical coordinates
                    fangles = @(phi1,phi2) f([ ...
                    r.*sin(phi1).*sin(phi2); ...
                    r.*cos(phi1).*sin(phi2); ...
                    r.*cos(phi2); ...
                    ]);

                    g = @(phi1,phi2) reshape(fangles(phi1(:)',phi2(:)').*sin(phi2(:)'),size(phi1)); % volume correcting term

                    dist = 0.5*integral2(g, 0, 2*pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                case 4
                    f = @(x) abs(this.pdf(x)-other.pdf(x));
                    r=1;

                    % hyperspherical coordinates
                    fangles = @(phi1,phi2,phi3) f([ ...
                    r.*sin(phi1).*sin(phi2).*sin(phi3); ...
                    r.*cos(phi1).*sin(phi2).*sin(phi3); ...
                    r.*cos(phi2).*sin(phi3); ...
                    r.*cos(phi3)
                    ]);

                    g = @(phi1,phi2,phi3) fangles(phi1,phi2,phi3) .* sin(phi2).*(sin(phi3)).^2; % volume correcting term
                    ga = @(phi1,phi2,phi3) reshape(g(phi1(:)', phi2(:)', phi3(:)'), size(phi1));

                    dist = 0.5*integral3(ga, 0, 2*pi, 0, pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                otherwise
                    error('Numerical calculation of total variation distance is currently not supported for this dimension.')
            end
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

