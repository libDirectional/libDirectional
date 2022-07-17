classdef (Abstract) AbstractHypercylindricalDistribution < AbstractLinPeriodicDistribution
    % Base class for distributions on the hypercylinder [0,2pi)^boundD x R^linD 
    
    methods
        function this = AbstractHypercylindricalDistribution()
            this.periodicManifoldType = 'hypertorus';
        end
        function p = plot(this, varargin)
            if this.boundD ~= 1
                error('not supported')
            end
            linSize = 3;
            switch this.linD
                case 1
                    % Creates a three-dimensional surface plot
                    step = 2*pi/100;
                    sigma = sqrt(this.linearCovariance());
                    m = this.mode();
                    [x,theta] = meshgrid(0:step:2*pi,linspace(m(1)-linSize*sigma,m(1)+linSize*sigma,100));
                    f = this.pdf([x(:)'; theta(:)']);
                    xlim([0,2*pi]);
                    f = reshape(f,size(x,1), size(x,2));
                    p = surf(x, theta, f, varargin{:});
                case 2
                    % Creates a three-dimensional plot
                    stepCirc = 0.5;
                    P = this.linearCovariance();
                    m = this.mode();
                    x1min = m(1)-linSize*sqrt(P(1,1));
                    x1max = m(1)+linSize*sqrt(P(1,1));
                    x2min = m(2)-linSize*sqrt(P(2,2));
                    x2max = m(2)+linSize*sqrt(P(2,2));
                    [X,Y,Z] = sphere(4);
                    clf 
                    hold on
                    color = jet;
                    [gridx,gridy,gridz]=ndgrid(0:stepCirc:2*pi, linspace(x1min,x1max,12), linspace(x2min,x2max,12));
                    fgrid=reshape(this.pdf([gridx(:)';gridy(:)';gridz(:)']),size(gridx));
                    fmax=max(fgrid(:));
                    sizes=0.5*stepCirc*fgrid/fmax;
                    arrayfun(@(x,y,z,currSize)surf(currSize*X+x,currSize*Y+y,currSize*Z+z, 'facecolor', color(1+floor(currSize*126/stepCirc),:)),...
                    gridx(sizes>0.01),gridy(sizes>0.01),gridz(sizes>0.01),sizes(sizes>0.01)); % Only use grid points for which size>0.01                    
                    hold off
                    xlabel('x_1')
                    ylabel('x_2')
                    zlabel('x_3')
                    setupAxisCircular('x')
                    view(40,20)
                    grid
                otherwise
                    error('not supported')
            end
        end
        
        function C = linearCovariance(this, approximateMean)
            % Overload this if non-numerical solutions exists
            arguments
                this (1,1) AbstractHypercylindricalDistribution
                approximateMean (:,1) double = NaN(this.linD,1)
            end
            assert(size(approximateMean,1)==this.linD);
            C = this.linearCovarianceNumerical(approximateMean);
        end
        
        function C = linearCovarianceNumerical(this, approximateMean)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
                approximateMean (:,1) double = NaN(this.linD,1)
            end
            assert(size(approximateMean,1)==this.linD);
            if anynan(approximateMean)
                approximateMean = this.linearMean();
            end
            if this.boundD==1 && this.linD==1
                C = integral2(@(x,y)(y-approximateMean).^2.*reshape(this.pdf([x(:)';y(:)']),size(x)),0,2*pi,-inf,inf);
            elseif this.boundD==2 && this.linD==1
                C = integral3(@(x,y,z)(z-approximateMean).^2.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,-inf,inf);
            elseif this.boundD==1 && this.linD==2
                C = NaN(2,2);
                C(1,1) = integral3(@(x,y,z)(y-approximateMean(1)).^2.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,-inf,inf,-inf,inf);
                C(1,2) = integral3(@(x,y,z)(y-approximateMean(1)).*(z-approximateMean(2)).*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,-inf,inf,-inf,inf);
                C(2,1) = C(1,2);
                C(2,2) = integral3(@(x,y,z)(z-approximateMean(2)).^2.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,-inf,inf,-inf,inf);
            else
                error('Cannot determine linear covariance for this dimension.');
            end
        end
        
        function mu = linearMean(this)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
            end
            mu = this.linearMeanNumerical;
        end
        
        function mu = linearMeanNumerical(this)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
            end
            if this.linD==1 && this.boundD==1
                % Only calculate for linear dimension
                mu = integral2(@(x,y)y.*reshape(this.pdf([x(:)';y(:)']),size(x)),0,2*pi,-inf,inf);
            elseif this.boundD==2 && this.linD==1
                mu = integral3(@(x,y,z)z.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,0,2*pi,-inf,inf);
            elseif this.boundD==1 && this.linD==2
                mu = NaN(2,1);
                mu(1) = integral3(@(x,y,z)y.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,-inf,inf,-inf,inf);
                mu(2) = integral3(@(x,y,z)z.*reshape(this.pdf([x(:)';y(:)';z(:)']),size(x)),0,2*pi,-inf,inf,-inf,inf);
            else
                error('Cannot determine linear covariance for this dimension.');
            end
        end
        
        function h = plotCylinder(this, limitsLinear)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
                limitsLinear (1,2) double = NaN(1,2)
            end
            assert(this.boundD == 1 && this.linD == 1);
            
            if anynan(limitsLinear)
                scaleLin = 3;
                m = this.linearMean();
                if ~anynan(m)
                    P = this.linearCovariance();
                end
                if ~anynan(m) && ~anynan(P)
                    limitsLinear(1) = m(1)-scaleLin*sqrt(P(1,1));
                    limitsLinear(2) = m(1)+scaleLin*sqrt(P(1,1));
                else
                    % Okay, we give up on this, just sample to find a
                    % suitable range.
                    s = this.sample(100);
                    limitsLinear=minmax(s(this.boundD+1:end));
                end
            end
            
            phi = linspace(0, 2*pi, 100);
            l = linspace(limitsLinear(1), limitsLinear(2), 100);
            [Phi,L] = meshgrid(phi,l);
            C = this.pdf([Phi(:)'; L(:)']);
            h = surf(cos(Phi),sin(Phi), L, reshape(C, size(Phi)));
            shading interp
        end
        
        function result = integral(this, l, r, replaceInfIfUsingNumerical)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
                l (:,1) double = [zeros(this.boundD,1); -inf(this.linD,1)]
                r (:,1) double = [2*pi*ones(this.boundD,1); inf(this.linD,1)]
                replaceInfIfUsingNumerical (1,1) double = false % Because it is expensive to obtain the limits
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
                this (1,1) AbstractHypercylindricalDistribution
                l (:,1) double = [zeros(this.boundD,1); -inf(this.linD,1)]
                r (:,1) double = [2*piones(this.boundD,1); inf(this.linD,1)]
                replaceInf (1,1) double = false % Because it is expensive to obtain the limits
            end
            assert(size(l,1)==this.dim);
            assert(size(r,1)==this.dim);    
            if replaceInf && any(isinf(l))||any(isinf(r))
                [lNonInf,rNonInf] = this.getIntegrationLimits();
                l(isinf(l))=lNonInf(isinf(l));
                r(isinf(r))=rNonInf(isinf(r));
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
                    error('Numerical moment calculation for this dimension is currently not supported');
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
                this (1,1) AbstractHypercylindricalDistribution
                scalingFactor (1,1) double {mustBePositive} = 10
            end
            
            l = zeros(this.linD + this.boundD, 1);
            r = 2*pi*ones(this.linD + this.boundD, 1);
            P = this.linearCovariance();
            m = this.mode();
            for i=this.boundD+1:this.boundD+this.linD % Change for linear dimensions
                l(i) = m(i)-scalingFactor*sqrt(P(i-this.boundD,i-this.boundD));
                r(i) = m(i)+scalingFactor*sqrt(P(i-this.boundD,i-this.boundD));
            end
        end
        
        function val = getManifoldSize(this)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
            end
            if this.linD==0 % If no linear dimension: Is hypertorus
                val = HypertoroidalUniformDistribution(this.boundD).getManifoldSize;
            else
                warning('Hypercylinder is of infinite size.');
                val = inf;
            end
        end
        
        function dist = conditionOnLinear(this, input, normalize)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
                input (:,1) double {mustBeNonempty}
                normalize (1,1) logical = true
            end
            assert(size(input,1) == this.linD);
            fCondUnnorm = @(x)this.pdf([x;repmat(input,1,size(x,2))]);
            dist = CustomHypertoroidalDistribution(fCondUnnorm, this.boundD);
            if normalize % Conditional need not be normalized
                dist = dist.normalize();
            end
        end
        
        function dist = conditionOnPeriodic(this, input, normalize)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
                input (:,1) double {mustBeNonempty}
                normalize (1,1) logical = true
            end
            assert(size(input,1) == this.boundD);
            fCondUnnorm = @(x)this.pdf([repmat(input,1,size(x,2));x]);
            dist = CustomLinearDistribution(fCondUnnorm, this.linD);
            if normalize % Conditional need not be normalized
                dist = dist.normalize();
            end
        end
        
        function m = mode(this)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
            end
            m = this.modeNumerical();
        end
        
        function m = modeNumerical(this, startingPoint)
            arguments
                this (1,1) AbstractHypercylindricalDistribution
                startingPoint (:,1) double = [pi*ones(this.boundD,1);zeros(this.linD,1)] % Start in the middle of the interval given as bound
            end
            m = fmincon(@(x)-this.pdf(x),startingPoint,[],[],[],[],...
                [zeros(this.boundD,1);-inf(this.linD,1)],[2*pi*ones(this.boundD,1);inf(this.linD,1)],[]);
            if isequal(m,startingPoint)
                warning('ModeNumerical:StoppedEarly','Mode was at starting point. This may indicate the optimizer stopped early.')
            end
        end
        
        function mu = hybridMean(this)
            % Do not fall back to hybridMeanNumerical because it should - whenver
            % available - use analytical formulations of hybridMoment
            m = this.hybridMoment();
            mu = [mod(atan2(m(2:2:2*this.boundD),m(1:2:2*this.boundD)),2*pi);m(2*this.boundD+1:end)];
        end
        
        function mu = hybridMeanNumerical(this)
            m = this.hybridMomentNumerical();
            mu = [mod(atan2(m(2:2:2*this.boundD),m(1:2:2*this.boundD)),2*pi);m(2*this.boundD+1:end)];
        end
        
        function s = sampleMetropolisHastings(this, n, proposal, startPoint, burnin, skipping)
            arguments
                this (1,1) AbstractDistribution
                n (1,1) {mustBePositive,mustBeInteger}
                % Proposal is HypercylindricalWN with zero mean and
                % identity matrix. Because there is no correlation, using randn
                % can be used.
                proposal (1,1) function_handle = @(x) [mod(x(1:this.boundD)+randn(this.boundD,1),2*pi);x(this.boundD+1:end)+randn(this.linD,1)]
                startPoint (1,1) = this.hybridMean()
                burnin (1,1) double = 10
                skipping (1,1) double = 5
            end
            s = sampleMetropolisHastings@AbstractDistribution(this, n, proposal, startPoint, burnin, skipping);
        end
        
        function mu = hybridMoment(this)
            % Do not fall back to numerical because it should - whenver
            % possible - use analytical formulations of hybridMoment
            mu = this.hybridMomentNumerical();
        end
        
        function m = hybridMomentNumerical(this)
            % Calculates mean of [cos(x1), sin(x1), cos(x2), .., cos(x_boundD), sin(x_boundD),
            % x_(boundD+1, x_(boundD+2), ..., x_(linD+boundD)]
            %
            % Returns:
            %   mu (linD + 2*boundD)
            %       Calculates mean of [cos(x1), sin(x1), cos(x2), .., cos(x_boundD), sin(x_boundD),
            %       x_(boundD+1, x_(boundD+2), ..., x_(linD+boundD)]
            
            % The linear part should be integrated from -Inf to Inf, but
            % Matlab's numerical integration does not handle that well.
            % For this reason we integrate from mu-10*sigma to mu+10*sigma,
            % which contains almost the entire probability mass
                        
            if this.boundD ~= 1
                error('not supported')
            end
           
            switch this.linD
                case 1
                    [l,r] = this.getIntegrationLimits();
                    
                    f1 = @(x,y) reshape(cos(x(:)').*this.pdf([x(:)';y(:)']), size(x,1), size(x,2));
                    f2 = @(x,y) reshape(sin(x(:)').*this.pdf([x(:)';y(:)']), size(x,1), size(x,2));
                    f3 = @(x,y) reshape(y(:)'.*this.pdf([x(:)';y(:)']), size(x,1), size(x,2));
                    
                    m1 = integral2(f1, l(1), r(1), l(2), r(2));
                    m2 = integral2(f2, l(1), r(1), l(2), r(2));
                    m3 = integral2(f3, l(1), r(1), l(2), r(2));
                    m = [m1; m2; m3];
                case 2
                    % todo this is extremely slow
                    [l,r] = this.getIntegrationLimits();

                    f1 = @(x,y,z) reshape(cos(x(:)').*this.pdf([x(:)';y(:)';z(:)']), size(x,1), size(x,2));
                    f2 = @(x,y,z) reshape(sin(x(:)').*this.pdf([x(:)';y(:)';z(:)']), size(x,1), size(x,2));
                    f3 = @(x,y,z) reshape(y(:)'.*this.pdf([x(:)';y(:)';z(:)']), size(x,1), size(x,2));
                    f4 = @(x,y,z) reshape(z(:)'.*this.pdf([x(:)';y(:)';z(:)']), size(x,1), size(x,2));
                    
                    m1 = integral3(f1, l(1), r(1), l(2), r(2), l(3), r(3));
                    m2 = integral3(f2, l(1), r(1), l(2), r(2), l(3), r(3));
                    m3 = integral3(f3, l(1), r(1), l(2), r(2), l(3), r(3));
                    mu4 = integral3(f4, l(1), r(1), l(2), r(2), l(3), r(3));
                    m = [m1; m2; m3; mu4];
                otherwise
                    error('unnsported')
            end            
        end
    end
    
end

