classdef AbstractHypersphericalDistribution
    % Abstract base class for distributions on the hypershere (S^dim)
    
    properties
        dim       % Dimension  (dim=2 => circle, dim=3 => sphere, ...)
    end
    
    methods (Abstract)
        % Evaluate pdf at positions stored in xa
        pdf(this, xa);
    end
    
    methods
        function plot(this, faces, gridFaces)
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
                    plot(phi, p);
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
                    if gridFaces > 0
                        surf(xSphereOuter, ySphereOuter, zSphereOuter, max(max(cSphere))*ones(size(xSphereOuter)), 'FaceColor', 'none');
                        hold on;
                    end
                    surf(xSphereInner, ySphereInner, zSphereInner, cSphere,'EdgeColor', 'none');
                    axis equal
                    colorbar
                    if ~holdStatus
                        hold off
                    end
                otherwise
                    error('cannot plot hyperspherical distribution with this number of dimensions');
            end
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
                r = randn(this.dim,n); % normal distribution
                r = r./repmat(sqrt(sum(r.^2)),this.dim,1); % normalize on unit sphere
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
        
        function s = sample(this, n)
            % Stocahastic sampling
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
            proposal = @(x) normalize(x + mvnrnd(zeros(this.dim,1),eye(this.dim))'); 
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
            %       distribution to compare to
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
            surfaceArea = 2 * pi^((dimension)/2)/gamma((dimension)/2); 
        end
        
    end 
end

