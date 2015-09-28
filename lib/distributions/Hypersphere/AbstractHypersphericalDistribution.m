classdef AbstractHypersphericalDistribution
    % Abstract base class for distributions on the hypershere (Sd)
    
    properties
        d       % Dimension  (d=2 => circle, d=3 => sphere, ...)
    end
    
    methods (Abstract)
        % Evaluate pdf at positions stored in xa
        pdf(this, xa);
    end
    
    methods
        function plot(this, faces, gridFaces)
        % Plots pdf of hyperspherical Distribution
        %   Parameters:
        %       faces - Number of faces for 3D Plot (default 100).
        %       gridFaces - Number of grid faces for 3D Plot (default 20, 0 to disable).
        
            switch this.d
                case 2
                    % use polar coordinates, plot angle->pdf(angle)
                    phi = 0:0.02:2*pi;
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
            if this.d==2
                %use matlab integration
                f = @(phi) this.pdf([cos(phi); sin(phi)]);
                i = integral(f,0,2*pi, 'AbsTol', 0.01);
            elseif this.d ==3
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
            elseif this.d==4     
                % use matlab integration
                f = @(x) this.pdf(x);
                r=1;
                
                % hyperspherical coordinates
                fangles = @(phi1,phi2,phi3) f([ ...
                r*sin(phi1)*sin(phi2)*sin(phi3); ...
                r*cos(phi1)*sin(phi2)*sin(phi3); ...
                r*cos(phi2)*sin(phi3); ...
                r*cos(phi3)
                ]);

                g = @(phi1,phi2,phi3) fangles(phi1,phi2,phi3) * sin(phi2)*(sin(phi3))^2; % volume correcting term
                ga = @(phi1,phi2,phi3) arrayfun(g, phi1, phi2, phi3);

                i = integral3(ga, 0, 2*pi, 0, pi, 0, pi, 'AbsTol', 1e-3, 'RelTol', 1e-3);
            else
                % use monte carlo integration
                n = 10000; % number of samples for integration
                r = randn(this.d,n); % normal distribution
                r = r./repmat(sqrt(sum(r.^2)),this.d,1); % normalize on unit sphere
                p = this.pdf(r);
                Sd = AbstractHypersphericalDistribution.computeUnitSphereSurface(this.d);
                i = sum(p)/n * Sd; % average value of pdf times surface area of unit sphere
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

