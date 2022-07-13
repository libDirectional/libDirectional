classdef (Abstract) AbstractHypersphereSubsetDistribution < AbstractPeriodicDistribution
    % For any distribution that is a subset (including equality) of the
    % hypersphere.
    % Parent class for AbstractHypersphericalDistribution and 
    % AbstractHyperhemisphericalDistributionto allow them to share code.
    
    methods
        function mu = meanDirection(this)
            % Calculate mean direction of pdf
            % Returns:
            %   mu (vector)
            %       mean direction
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
            end
            mu = meanDirectionNumerical(this);
        end

        function mu = meanDirectionNumerical(this, integrationBoundaries)
            % Calculate mean direction of pdf. Specify part of hypersphere
            % by setting integrationBoundaries (this.dim x 2 matris)
            % accordingly
            % Returns:
            %   mu (vector)
            %       mean direction
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
                integrationBoundaries (:,2) double
            end
            mu = NaN(this.dim,1);
            switch this.dim
                case 2
                    for i=1:2
                        f = @(x) x(i,:).*this.pdf(x);
                        fAngles = @(phi) reshape(f([cos(phi(:)'); sin(phi(:)')]),size(phi));
                        mu(i) = integral(fAngles,integrationBoundaries(1,1), integrationBoundaries(1,2), 'AbsTol', 0.01);
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
                        mu(i) = integral2(g, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                            integrationBoundaries(2,1),integrationBoundaries(2,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                        mu(i) = integral3(ga, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                            integrationBoundaries(2,1),integrationBoundaries(2,2),...
                            integrationBoundaries(3,1),integrationBoundaries(3,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
                    end
                otherwise
                    error('Unsupported');
            end
            if norm(mu)<1e-9
                warning('Density may not have actually have a mean direction because integral yields a point very close to the origin.')
            end
            mu = mu/norm(mu);
        end

        function m = moment(this)
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
            end
            m = this.momentNumerical();
        end

        function m = momentNumerical(this, integrationBoundaries)
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
                integrationBoundaries (:,2) double
            end
            m = NaN(this.dim);
            fGen = @(i,j) @(points) reshape(this.pdf(points).*points(i,:).*points(j,:), [1,size(points,2)]);
            r = 1;
            switch this.dim
                case 2
                    for i=1:2
                        for j=1:2
                            fCurr = fGen(i,j);
                            fangles = @(phi) fCurr([sin(phi); cos(phi)]);
                            g = @(phi)reshape(fangles(phi(:)'),size(phi));
                            m(i,j) = integral(g,integrationBoundaries(1,1), integrationBoundaries(1,2));
                        end
                    end
                case 3
                    for i=1:3
                        for j=1:3
                            fCurr = fGen(i,j);
                            fangles = @(phi1,phi2) fCurr([ ...
                            r.*sin(phi1).*sin(phi2); ...
                            r.*cos(phi1).*sin(phi2); ...
                            r.*cos(phi2); ...
                            ]);
        
                            g = @(phi1,phi2) reshape(fangles(phi1(:)',phi2(:)').*sin(phi2(:)'),size(phi1)); % volume correcting term
                            m(i,j) = integral2(g, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                                integrationBoundaries(2,1), integrationBoundaries(2,2));
                        end
                    end
                case 4
                    for i=1:4
                        for j=1:4
                            fCurr = fGen(i,j);
                            fangles = @(phi1,phi2,phi3) fCurr([ ...
                                r.*sin(phi1).*sin(phi2).*sin(phi3); ...
                                r.*cos(phi1).*sin(phi2).*sin(phi3); ...
                                r.*cos(phi2).*sin(phi3); ...
                                r.*cos(phi3)
                            ]);
                            
                            g = @(phi1,phi2,phi3) fangles(phi1,phi2,phi3) .* sin(phi2).*(sin(phi3)).^2; % volume correcting term
                            ga = @(phi1,phi2,phi3) reshape(g(phi1(:)', phi2(:)', phi3(:)'), size(phi1));
                            
                            m(i,j) = integral3(ga, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                                integrationBoundaries(2,1), integrationBoundaries(2,2),...
                                integrationBoundaries(3,1), integrationBoundaries(3,2));
                        end
                    end
                otherwise
                    error('Dimension not supported.')
            end
        end

        function m = meanAxis(this)
            % Do not just fall back to .momentNumerical() per default
            % because the class may provide a good function to generate
            % .moment().
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
            end
            mom = this.moment();
            [V, D] = eig(mom);
            [Dsorted, order] = sort(diag(D));
            Vsorted = V(:,order);
            if abs(Dsorted(end)/Dsorted(end-1))<1.01
                warning('Eigenvalues are very similar. Axis may be unreliable.');
            end            
            if Vsorted(end,end)>=0
                m = Vsorted(:,end);
            else
                m = -Vsorted(:,end);
            end
        end

        function m = meanAxisNumerical(this)
            % This function calls .momentNumerical() for numerical
            % computation of the moment.
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
            end
            mom = this.momentNumerical();
            [V, D] = eig(mom);
            [Dsorted, order] = sort(diag(D));
            Vsorted = V(:,order);
            if abs(Dsorted(end)/Dsorted(end-1))<1.01
                warning('Eigenvalues are very similar. Axis may be unreliable.');
            end            
            if Vsorted(end,end)>=0
                m = Vsorted(:,end);
            else
                m = -Vsorted(:,end);
            end
        end

        function i = integral(this)
            % Calculate integral of pdf to check normalization
            % (should always be 1)
            % Returns:
            %   i (scalar)
            %       integral over hypersphere surface of pdf (uses
            %       approximation, not very accurate for higher dimensions)
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
            end
            i = this.integralNumerical();
        end

        function i = integralNumerical(this, integrationBoundaries)
            % Calculate integral of pdf to check normalization
            % (should always be 1)
            % Returns:
            %   i (scalar)
            %       integral over hypersphere surface of pdf (uses
            %       approximation, not very accurate for higher dimensions)
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
                integrationBoundaries (:,2) double
            end
            switch this.dim
                case 2
                    %use matlab integration
                    f = @(phi) this.pdf([cos(phi); sin(phi)]);
                    i = integral(f,integrationBoundaries(1,1), integrationBoundaries(1,2), 'AbsTol', 0.01);
                case 3
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
    
                    i = integral2(g, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                        integrationBoundaries(2,1), integrationBoundaries(2,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
                case 4   
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
    
                    i = integral3(ga, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                        integrationBoundaries(2,1), integrationBoundaries(2,2),...
                        integrationBoundaries(3,1), integrationBoundaries(3,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
            otherwise
                    error('Dimension not supported.')
            end
        end

        function result = entropy(this)
            % Calculates the entropy analytically if possible, 
            % fall back to numerical calculation by default
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
            end
            result = this.entropyNumerical();
        end

        function i = entropyNumerical(this, integrationBoundaries)
            % Calculates the entropy numerically
            %
            % Returns:
            %   i (scalar)
            %       entropy of the distribution
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
                integrationBoundaries (:,2) double
            end
            if this.dim==2
                %use matlab integration
                f = @(phi) this.pdf([cos(phi); sin(phi)]).*log(this.pdf([cos(phi); sin(phi)]));
                i = -integral(f,integrationBoundaries(1,1), integrationBoundaries(1,2), 'AbsTol', 0.01);
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

                i = -integral2(ga, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                        integrationBoundaries(2,1), integrationBoundaries(2,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                i = -integral3(ga, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                            integrationBoundaries(2,1),integrationBoundaries(2,2),...
                            integrationBoundaries(3,1),integrationBoundaries(3,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
            else
                error('not supported')
            end
        end

        function m = mode(this)
            % Calculate the mode by finding maximum of PDF numerically 
            % Returns:
            %   m (column vector)
            %       mode of the distribution 
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
            end
            m = this.modeNumerical(); % fallback to numerical calculation (nonlinear optimization)
        end

        function dist = hellingerDistanceNumerical(this, other, integrationBoundaries)
            % Numerically calculates the Hellinger distance to another
            % distribution.
            %
            % Parameters:
            %   other (AbstractHypersphericalDistribution)
            %       distribution to compare to
            % Returns:
            %   dist (scalar)
            %       hellinger distance of this distribution to other distribution
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
                other (1,1) AbstractHypersphereSubsetDistribution
                integrationBoundaries (:,2) double
            end
            assert(this.dim==other.dim,'Cannot compare distributions with different number of dimensions');
            
            % Implementation is always performed using matlab integration
            % (Implementation is similar to .integral)
            switch this.dim
                case 2
                    f = @(phi) (sqrt(this.pdf([cos(phi); sin(phi)]))-sqrt(other.pdf([cos(phi); sin(phi)]))).^2;
                    dist = 0.5*integral(f,integrationBoundaries(1,1), integrationBoundaries(1,2), 'AbsTol', 0.01);
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

                    dist = 0.5*integral2(g, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                        integrationBoundaries(2,1), integrationBoundaries(2,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                    dist = 0.5*integral3(ga, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                            integrationBoundaries(2,1),integrationBoundaries(2,2),...
                            integrationBoundaries(3,1),integrationBoundaries(3,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
                otherwise
                    error('Numerical calculation of Hellinger distance is currently not supported for this dimension.')
            end
        end
        
        function dist = totalVariationDistanceNumerical(this, other, integrationBoundaries)
            % Numerically calculates the total varation distance to another distribution
            %
            % Parameters:
            %   other (AbstractHypersphericalDistribution)
            %       distribution to compare with
            % Returns:
            %   dist (scalar)
            %       total variation distance of this distribution to other distribution
            arguments
                this (1,1) AbstractHypersphereSubsetDistribution
                other (1,1) AbstractHypersphereSubsetDistribution
                integrationBoundaries (:,2) double
            end
            assert(this.dim==other.dim, 'Cannot compare distributions with different number of dimensions');
            
            switch this.dim
                case 2
                    f = @(phi) abs(this.pdf([cos(phi); sin(phi)])-other.pdf([cos(phi); sin(phi)]));
                    dist = 0.5*integral(f,integrationBoundaries(1,1), integrationBoundaries(1,2), 'AbsTol', 0.01);
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

                    dist = 0.5*integral2(g, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                            integrationBoundaries(2,1),integrationBoundaries(2,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
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

                    dist = 0.5*integral3(ga, integrationBoundaries(1,1), integrationBoundaries(1,2),...
                            integrationBoundaries(2,1),integrationBoundaries(2,2),...
                            integrationBoundaries(3,1),integrationBoundaries(3,2), 'AbsTol', 1e-3, 'RelTol', 1e-3);
                otherwise
                    error('Numerical calculation of total variation distance is currently not supported for this dimension.')
            end
        end
    end 
end

