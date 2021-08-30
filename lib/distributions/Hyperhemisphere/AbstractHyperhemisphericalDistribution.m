classdef (Abstract) AbstractHyperhemisphericalDistribution < AbstractDistribution
    % Abstract base class for distributions on the hypershere (S^dim)
    
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
                    phi = linspace(0,pi,320);
                    x = [cos(phi); sin(phi)];
                    p = this.pdf(x);
                    h = plot(phi, p);
                case 3
                    % plot sphere, pdf is represented by color on sphere
                    
                    if nargin < 2
                        faces = 100;
                    else
                        assert(mod(faces,2)==0);
                    end
                    
                    if nargin < 3
                        gridFaces = 20;
                    else
                        assert(mod(gridFaces,2)==0);
                    end
                    
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
                    error('cannot plot hyperspherical distribution with this number of dimensions');
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
            warning(['The result is the mean direction on the upper hemisphere ',...
            'along the last dimension. It is not a mean of a symmetric ',...
            'distribution, which would not have a proper mean. ',...
            'It is also not one of the modes of the symmetric distribution ',...
            'since it is biased toward [0;...;0;1] because the lower half ',...
            'is considered inexistant.']);
            mu = NaN(this.dim,1);
            switch this.dim
                case 2
                    for i=1:2
                        f = @(x) x(i,:).*this.pdf(x);
                        fAngles = @(phi) reshape(f([cos(phi(:)'); sin(phi(:)')]),size(phi));
                        mu(i) = integral(fAngles,0,pi, 'AbsTol', 0.01);
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
                        mu(i) = integral2(g, 0, 2*pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                        mu(i) = integral3(ga, 0, 2*pi, 0, pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                    end
                otherwise
                    % use monte carlo integration
                    Sd = 0.5*AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim);
                    
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
                i = integral(f,0,pi, 'AbsTol', 0.01);
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

                i = integral2(g, 0, 2*pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                i = integral3(ga, 0, 2*pi, 0, pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
            else
                % use monte carlo integration
                n = 10000; % number of samples for integration
                r = HyperhemisphericalUniformDistribution.sample(n);
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
                i = -integral(f,0,pi, 'AbsTol', 0.01);
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

                i = -integral2(ga, 0, 2*pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                i = -integral3(ga, 0, 2*pi, 0, pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
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
                    dist = 0.5*integral(f,0,pi, 'AbsTol', 0.01);
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

                    dist = 0.5*integral2(g, 0, 2*pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                    dist = 0.5*integral3(ga, 0, 2*pi, 0, pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                otherwise
                    error('Numerical calculation of Hellinger distance is currently not supported for this dimension.')
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
                    dist = 0.5*integral(f,0,pi, 'AbsTol', 0.01);
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

                    dist = 0.5*integral2(g, 0, 2*pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                    dist = 0.5*integral3(ga, 0, 2*pi, 0, pi, 0, pi/2, 'AbsTol', 1e-3, 'RelTol', 1e-3);
                otherwise
                    error('Numerical calculation of total variation distance is currently not supported for this dimension.')
            end
        end
        function s = getManifoldSize(this)
            s = 0.5 * AbstractHypersphericalDistribution.computeUnitSphereSurface(this.dim);
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

