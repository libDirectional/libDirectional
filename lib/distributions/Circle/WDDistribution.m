classdef WDDistribution < AbstractCircularDistribution
    % Wrapped Dirac distribution (i.e., discrete distribution)
    % with Dirac positions d and weights w
    %
    % see Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Recursive Nonlinear Filtering for Angular Data Based on Circular Distributions
    % Proceedings of the 2013 American Control Conference (ACC 2013), Washington D. C., USA, June 2013
    
    properties
        d
        w
    end
    
    methods
        function this = WDDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (1 x L)
            %       Dirac locations in [0,2pi)
            %   w_ (1 x L)
            %       weights for each Dirac            
            assert(all(size(d_,1)==1), 'd must be a row vector');
            this.d = mod(d_,2*pi);
            if (nargin<2)
                %all diracs have equal weights
                this.w = ones(size(this.d))/length(this.d);
            else
                assert(all(size(d_)==size(w_)), 'sizes of d_ and w_ must match');
                assert(all(w_ >=0), 'weights must be non-negative');
                assert(sum(w_)>0, 'sum of weights must be larger than zero');
                this.w = w_/sum(w_);
            end
        end
        
        function p = pdf(this, xa)
            % Placeholder, pdf does not exist for wrapped Dirac distributions
            p = 0; %WD does not have a proper pdf
            warning('PDF:UNDEFINED', 'pdf is not defined')
        end
        
        function val=cdf(this,xa,startingPoint)
            % Evaluate cumulative distribution function 
            %
            % Parameters:
            %   xa (1 x n)
            %       points where the cdf should be evaluated
            %   startingPoint (scalar)
            %       point where the cdf is zero (starting point can be
            %       [0,2pi) on the circle, default 0
            % Returns:
            %   val (1 x n)
            %       cdf evaluated at columns of xa
            if nargin <=2
                startingPoint = 0;
            end
            % shift d and xa so that startingPoint is beginning of the 2pi
            % interval
            dshifted=this.d+(this.d<startingPoint)*2*pi;
            xamod=mod(xa,2*pi);
            xamodshift=xamod+(xamod<startingPoint)*2*pi; 
            val=arrayfun(@(xi)sum(this.w(dshifted<=xi)),xamodshift);
        end
                
        function m = trigonometricMoment(this,n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            m = sum(exp(1i*n*this.d).*this.w);
        end
        
        function wd=shift(this,angle)
            wd=WDDistribution(this.d+angle,this.w); % Mod is done in constructor
        end
        
        function m = trigonometricMomentNumerical(this,n)
            % Disable numerical calculation of trigonometric moments since it relies on the pdf
            error('not supported');
        end
        
        function d = squaredDistanceNumerical(this, other)
            error('not supported');
        end
        
        function kld = kldNumerical(this, other)
            error('not supported');
        end
        
        function v = kuiperTest(this, other)
            % Evaluate Kuiper's test
            %
            % Parameters:
            %   other (AbstractCircularDistributuon)
            %       distribution to compare to
            % Returns:
            %   v (scalar)
            %       test result
            %
            % see
            % http://en.wikipedia.org/wiki/Kuiper's_test
            % and
            % Watson, G. S. Goodness-Of-Fit Tests on a Circle Biometrika, 
            % Biometrika Trust, 1961, 48, 109-114
            
            assert(isa(other, 'AbstractCircularDistribution'));
            
            [dsorted, order] = sort(this.d);
            wsorted = this.w(order);
            z = other.cdf(dsorted);
            dplus = max(cumsum(wsorted) - z);
            dminus = max(z - [0, cumsum(wsorted(1:end-1))]);
            v = dplus + dminus;
        end
                        
        function wn = toWNunwrappingEM(this)
            % Obtain a WN by unwrapping and EM
            %
            % Returns:
            %   wn (WNDistribution)
            %       wn distribution obtained by Fisher's algorithm
            %
            % "Time Series Analysis of Circular Data", N. I. Fisher and A.
            % J. Lee, Journal of the Royal Statistical Society. Series B (Methodological), 
            % Vol. 56, No. 2(1994), pp. 327-339
            kmax = 3;
            
            % start value
            mu=0;
            sigma=1;
            
            epsilon = 1E-7; % stopping condition
            maxIter = 30; % maximum number of iterations
            
            for iteration=1:maxIter
                samples=zeros(length(this.d),2*kmax+1);
                p=zeros(length(this.d),2*kmax+1);
                
                % E-Step
                for i=1:length(this.d)
                    for k=-kmax:kmax
                        % probability that sample i was wrapped k times
                        samples(i,k+kmax+1) = this.d(i)+2*k*pi;
                        p(i,k+kmax+1) = this.w(i)*normpdf(samples(i,k+kmax+1), mu, sigma);
                    end
                    % normalize per sample
                    p(i,:) =  p(i,:)/sum(p(i,:));                
                end

                % M-Step
                p = p(:)/sum(p(:)); % normalize all
                muOld = mu;
                sigmaOld = sigma;
                mu = sum(p(:).*samples(:));
                sigma = sqrt(sum(p(:).*(samples(:)-mu).^2));
                
                % check stopping condition
                if abs(mu - muOld) < epsilon && abs(sigma - sigmaOld) < epsilon
                    break
                end
            end
            
            wn = WNDistribution(mu,sigma);
        end
                
        function wd = applyFunction(this,f)
            % Apply a function f(x) to each Dirac component and obtain its new position
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi) to [0,2*pi)
            % Returns:
            %   wd (WDDistribution)
            %       distribution with new Dirac locations (and same
            %       weights as before)
            d_ = arrayfun(f, this.d);
            wd = WDDistribution(d_, this.w);
        end
        
        function wd = reweigh(this, f)
            % Uses a function f(x) to calculate the weight of each Dirac
            % component. The new weight is given by the product of the old 
            % weight and the weight obtained with f. Restores normalization
            % afterwards.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi) to [0, infinity)
            % Returns:
            %   wd (WDDistribution)
            %       distribution with new weights and same Dirac locations
            wNew = zeros(1, length(this.d));
            for i=1:length(this.d)
                wNew(i) = f(this.d(i));
            end
            wd = WDDistribution(this.d, wNew .* this.w);
        end
        
        function result = integralNumerical(this, l, r)
            error('not supported');
        end
        
        function result = integral(this, l, r)
            % Calculates the integral of the pdf from l to r analytically
            %
            % Parameters:
            %   l (scalar)
            %       left bound of integral, default 0
            %   r (scalar)
            %       right bound of integral, default 2*pi
            % Returns:
            %   result (scalar)
            %       value of the integral
            if nargin < 2
                l = 0;
            end
            if nargin < 3
                r = 2*pi;
            end
            
            if l<=r
                result = floor((r-l)/(2*pi)); %each full 2pi interval contributes 1 to the result because of normalization
                r = mod(r,2*pi);
                l = mod(l,2*pi);
                if l<= r
                    result = result + sum(this.w(this.d>=l & this.d<=r));
                else
                    result = result + 1 - sum(this.w(this.d>=r & this.d<=l));
                end
            else
                result = -this.integral(r,l);
            end
        end
        
        function h = plot2d(this, varargin)
            % Create a 2D plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot command
            % Returns:
            %   p (scalar)
            %       plot handle
            h = stem(this.d, this.w, varargin{:});
        end
        
        function h = plot2dcircular(this, varargin)
            % Create a circular 2D plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot command
            % Returns:
            %   p (scalar)
            %       plot handle
            error('not implemented');
        end
        
        function h = plot3d(this, varargin)
            % Create a 3D plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot command
            % Returns:
            %   p (1 x 2)
            %       plot handles for plot3 and stem3
            h = stem3(cos(this.d), sin(this.d), this.w, varargin{:});
        end
        
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (1 x n)
            %       n samples on the circle
            ids = discretesample(this.w,n);
            s = this.d(ids);
        end
        
        function s = sampleCdf(this, n)
            % Disable sampling algorithm relying on cdf
            error('not supported');
        end
        
        function s = sampleMetropolisHastings(this, n)
            % Disable sampling algorithm relying on pdf
            error('not supported');
        end
        
        function wd = convolve(this, wd2)
            % Calculates the convolution of this and wd2
            %
            % Parameters:
            %   wd2 (WDDistribution)
            %       distribution to convolve with
            % Returns:
            %   wd (WDDistribution)
            %       convolution of this and wd2
            assert(isa(wd2, 'WDDistribution'));
            
            L = length(this.d);
            Lw = length(wd2.d);
            d_ = zeros(1,L*Lw);
            w_ = zeros(1,L*Lw);
            for j=1:L
                for l=1:Lw;
                    w_(1,j+L*(l-1)) = this.w(j) * wd2.w(l);
                    d_(1,j+L*(l-1)) = this.d(j) + wd2.d(l);
                end
            end
            wd = WDDistribution(d_,w_);
        end
        
        function result = entropyNumerical(this)
            warning('entropy is not defined in a continous sense')
            result = 0;
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            warning('ENTROPY:DISCRETE','entropy is not defined in a continous sense')
            % return discrete entropy
            % The entropy only depends on the weights!
            result = -sum(this.w.*log(this.w));
        end
        
        function l = logLikelihood(this,samples)
            % Calculates the log-likelihood of the given samples
            %
            % Parameters:
            %   sampls (1 x n)
            %       n samples on the circle
            % Returns:
            %   l (scalar)
            %       log-likelihood of obtaining the given samples
            assert(size(samples,1)>=1);
            assert(size(samples,2)==1);
            
            error('unsupported')
        end
    end
    
end

