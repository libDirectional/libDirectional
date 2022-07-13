classdef (Abstract) AbstractLinearDistribution < AbstractDistribution
    methods
        function mu = mean(this)
            arguments
                this (1,1) AbstractLinearDistribution
            end
            mu = this.meanNumerical();
        end
        
        function mu = meanNumerical(this)
            arguments
                this (1,1) AbstractLinearDistribution
            end
            if this.dim==1
                mu = integral(@(x)x.*this.pdf(x),-inf,inf);
            elseif this.dim==2
                mu = NaN(2,1);
                mu(1) = integral2(@(x,y)x.*reshape(this.pdf([x(:)';y(:)']),size(x)),-inf,inf,-inf,inf);
                mu(2) = integral2(@(x,y)y.*reshape(this.pdf([x(:)';y(:)']),size(x)),-inf,inf,-inf,inf);
            elseif this.dim==3
                mu = NaN(3,1);
                mu(1) = integral3(@(x,y,z)x.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                mu(2) = integral3(@(x,y,z)y.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                mu(3) = integral3(@(x,y,z)z.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
            else
                error('Dimension currently not supported for all types of densities.');
            end
        end

        function C = covariance(this)
            arguments
                this (1,1) AbstractLinearDistribution
            end
            C = this.covarianceNumerical();
        end
        
        function C = covarianceNumerical(this)
            arguments
                this (1,1) AbstractLinearDistribution
            end
            mu = this.mean();
            if this.dim==1
                C = integral(@(x)(x-mu).^2.*this.pdf(x),-inf,inf);
             elseif this.dim==2
                C = NaN(2,2);
                C(1,1) = integral2(@(x,y)(x-mu(1)).^2.*reshape(this.pdf([x(:)';y(:)']),size(x)),-inf,inf,-inf,inf);
                C(1,2) = integral2(@(x,y)(x-mu(1)).*(y-mu(2)).*reshape(this.pdf([x(:)';y(:)']),size(x)),-inf,inf,-inf,inf);
                C(2,2) = integral2(@(x,y)(y-mu(2)).^2.*reshape(this.pdf([x(:)';y(:)']),size(x)),-inf,inf,-inf,inf);
                C(2,1) = C(1,2);
            elseif this.dim==3
                C = NaN(3,3);
                C(1,1) = integral3(@(x,y,z)(x-mu(1)).^2.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                C(1,2) = integral3(@(x,y,z)(x-mu(1)).*(y-mu(2)).*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                C(1,3) = integral3(@(x,y,z)(x-mu(1)).*(z-mu(3)).*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                C(2,2) = integral3(@(x,y,z)(y-mu(2)).^2.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                C(2,3) = integral3(@(x,y,z)(y-mu(2)).*(z-mu(3)).*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                C(3,3) = integral3(@(x,y,z)(z-mu(3)).^2.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),-inf,inf,-inf,inf,-inf,inf);
                C(2,1) = C(1,2);
                C(3,1) = C(1,3);
                C(3,2) = C(2,3);
            else
                error('Dimension currently not supported for all types of densities.');
            end
        end
        
        function m = mode(this)
            m = modeNumerical(this);
        end
        
        function m = modeNumerical(this, startingPoint)
            arguments
                this (1,1) AbstractLinearDistribution
                startingPoint (:,1) double = zeros(this.dim,1)
            end
            m = fmincon(@(x)-this.pdf(x),startingPoint,[],[],[],[],-inf(this.dim,1),inf(this.dim,1),[]);
            if isequal(m,startingPoint)
                warning('ModeNumerical:StoppedEarly','Mode was at starting point. This may indicate the optimizer stopped early.')
            end
        end
        
        function h = plot(this,range,varargin)
            arguments
                this (1,1) AbstractLinearDistribution
                range {mustBeNumeric} = []
            end
            arguments (Repeating)
                varargin
            end
            % We have to choose some heuristic regarding the area to plot.
            % This works pretty good for Gaussians.
            assert(isempty(range) || numel(range) == 2*this.dim);
            mu = this.mean();
            C = this.covariance();
            switch this.dim
                case 1
                    if isempty(range)
                        scaling = sqrt(chi2inv(0.99,this.dim));
                        range = [mu-scaling*C,mu+scaling*C];
                    end
                    h = fplot(@(x)this.pdf(x),range,varargin{:});
                case 2
                    if isempty(range)
                        scaling = sqrt(chi2inv(0.99,this.dim));
                        bounds = scaling*vecnorm(C);
                        range = [mu(1)-bounds(1),mu(1)+bounds(1),...
                            mu(2)-bounds(2), mu(2)+bounds(2)];
                    end
                    fun = @(x,y)reshape(this.pdf([x(:)';y(:)']),size(x));
                    h = fsurf(fun,range,varargin{:});
                otherwise
                    error('Dimension not supported');
            end
        end
        
        function val = getManifoldSize(this)
            arguments
                this (1,1) AbstractLinearDistribution
            end
            warning('R^%d is of infinite size.', this.dim);
            val = inf;
        end
        
        function s = sampleMetropolisHastings(this, n, proposal, startPoint, burnIn, skipping)
            arguments
                this (1,1) AbstractDistribution
                n (1,1) {mustBePositive,mustBeInteger}
                proposal (1,1) function_handle = @(x) x + randn(this.dim,1)
                startPoint (:,1) double = this.mean()
                burnIn (1,1) double = 10
                skipping (1,1) double = 5
            end
            s = sampleMetropolisHastings@AbstractDistribution(this, n, proposal, startPoint, burnIn, skipping);
        end
        
        function result = integral(this, l, r, replaceInfIfUsingNumerical)
            arguments
                this (1,1) AbstractLinearDistribution
                l (:,1) double = -inf(this.dim,1)
                r (:,1) double = inf(this.dim,1)
                replaceInfIfUsingNumerical (1,1) double = false
            end
            % Integrates the density to check normalization
            %
            % Returns:
            %   r (scalar)
            %       integral over the entire density
            result = this.integralNumerical(l,r, replaceInfIfUsingNumerical);
        end
        
        function result = integralNumerical(this, l, r, replaceInf)
            % Numerically calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l (dim x 1 column vector)
            %       left bound of integral in each dimension, default -inf
            %       or 0
            %   r (dim x 1 column vector)
            %       right bound of integral in each dimension, default 2*pi
            %       or inf
            arguments
                this (1,1) AbstractLinearDistribution
                l (:,1) double = -inf(this.dim,1)
                r (:,1) double = inf(this.dim,1)
                replaceInf (1,1) double = false
            end
            assert(size(l,1)==this.dim);
            assert(size(r,1)==this.dim);    
            if replaceInf && (any(isinf(l))||any(isinf(r)))
                [lNonInf,rNonInf] = this.getIntegrationLimits();
                l(isinf(l))=lNonInf(isinf(l));
                r(isinf(l))=rNonInf(isinf(r));
            end
            
            switch this.dim
                case 1
                    result = integral(@(x) this.pdf(x), l, r);
                case 2
                    f = @(x,y) reshape(this.pdf([x(:)';y(:)']), size(x));
                    result = integral2(f, l(1), r(1), l(2), r(2));
                case 3
                    f = @(x,y,z) reshape(this.pdf([x(:)';y(:)';z(:)']), size(x));
                    result = integral3(f, l(1), r(1), l(2), r(2), l(3), r(3));
                otherwise
                    error('Numerical integral calculation for this dimension is currently not supported');
            end
        end
        
        function [l,r] = getIntegrationLimits(this, scalingFactor)
            % Returns suggested limits for integration over the whole
            % density
            %            
            % The linear part should be integrated from -Inf to Inf, but
            % Matlab's numerical integration does not handle that well.
            % When we can obtain the covariance of the linear part easily,
            % we integrate from mu-10*sigma to mu+10*sigma,
            % which contains almost the entire probability mass. The
            % circular part is integrated form 0 to 2pi.
            %
            % Returns:
            %   l (linD+boundD x 1 column vector)
            %       lower integration bound
            %   r (linD+boundD x 1 column vector)
            %       upper integration bound
            arguments
                this (1,1) AbstractLinearDistribution
                scalingFactor (1,1) double {mustBePositive} = 10
            end
            P = this.linearCovariance();
            m = this.mode();
            l = NaN(this.dim,1);
            r = NaN(this.dim,1);
            for i=this.boundD+1:this.boundD+this.linD % Change for linear dimensions
                l(i) = m(i)-scalingFactor*sqrt(P(i-this.boundD,i-this.boundD));
                r(i) = m(i)+scalingFactor*sqrt(P(i-this.boundD,i-this.boundD));
            end
        end

        function lineHandles = plotState(this)
            arguments
                this (1,1) AbstractLinearDistribution
            end
            assert(this.dim==3);
            [x,y,z]=sphere(150); % Create smooth sphere
            h=mesh(x,y,z);
            

            skipx=10;
            skipy=10;
            x=get(h,'xdata'); % Get lines from smooth sphere
            y=get(h,'ydata');
            z=get(h,'zdata');
            delete(h)

            [V,D] = eig(this.C);
            allCoords=V*sqrt(D)*[x(:)';y(:)';z(:)'] + this.mean();
            x = reshape(allCoords(1,:),size(x));
            y = reshape(allCoords(2,:),size(y));
            z = reshape(allCoords(3,:),size(z));

            xKeep=x(1:skipx:end,:); % Only plot some of the lines as grid would be too fine otherwise
            yKeep=y(1:skipx:end,:);
            zKeep=z(1:skipx:end,:);
            lineHandles=line(xKeep',yKeep',zKeep','color',0.7*[1 1 1]);

            xKeep=x(:,1:skipy:end);
            yKeep=y(:,1:skipy:end);
            zKeep=z(:,1:skipy:end);
            lineHandles=[lineHandles;...
                line(xKeep,yKeep,zKeep,'color',0.7*[1 1 1])]';
        end
    end
end

