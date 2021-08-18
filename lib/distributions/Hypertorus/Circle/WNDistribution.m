classdef WNDistribution < AbstractCircularDistribution
    % wrapped normal distribution
    %
    % see Sreenivasa Rao Jammalamadaka and A. SenGupta, "Topics in Circular
    % Statistics", 2001, Sec. 2.2.6, page 44f.
    
    properties
        mu (1,1) double
        sigma (1,1) double {mustBeNonnegative}
    end
    
    methods
        function this = WNDistribution(mu_, sigma_)
            arguments
                mu_ (1,1) double
                sigma_ (1,1) double {mustBeNonnegative}
            end
            % Constructor
            this.mu = mod(mu_,2*pi);
            this.sigma = sigma_;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (1 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            arguments
                this (1,1) WNDistribution
                xa (1,:) double
            end
            p = wnpdf(xa, this.mu, this.sigma);
        end
        
        function val = cdf(this, xa, startingPoint, n)
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
            arguments
                this (1,1) WNDistribution
                xa (1,:) double
                startingPoint (1,1) double = 0
                n (1,1) double {mustBeInteger,mustBePositive} = 10  %Use n for predefined number of runs to use vectorization
            end
            startingPoint = mod(startingPoint,2*pi);
            xa = mod(xa,2*pi);
            
            ncdf = @(from,to) 1/2*(erf((this.mu-from)/(sqrt(2)*this.sigma))-erf((this.mu-to)/(sqrt(2)*this.sigma)));
            val = ncdf(startingPoint,xa);
            for i=1:n
                val = val+ncdf(startingPoint+2*pi*i,xa+2*pi*i)+ncdf(startingPoint-2*pi*i,xa-2*pi*i);
            end
            val = (xa<startingPoint)+val; %1+int(a,b) if a>b, otherweise int(a,b)
        end
        
        function m = trigonometricMoment(this, n)
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            arguments
                this (1,1) WNDistribution
                n (1,1) double {mustBeInteger}
            end
            
            m = exp(1i*n*this.mu - n^2*this.sigma^2/2);
        end
        
        function wn = multiplyVM(this,wn2)
            % Multiply two WN distributions (approximate by means of intermediate VM)
            %
            % Parameters:
            %   wn2 (WNDistribution)
            %       distribution to multiply with
            % Returns:
            %   wn (WNDistribution)
            %       product of this and wn2
            %
            % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
            % Recursive Nonlinear Filtering for Angular Data Based on Circular Distributions
            % Proceedings of the 2013 American Control Conference (ACC 2013), Washington D. C., USA, June 2013.
            arguments
                this (1,1) WNDistribution
                wn2 (1,1) AbstractCircularDistribution
            end
            
            vm1 = this.toVM();
            vm2 = wn2.toVM();
            vm = vm1.multiply(vm2);
            wn = vm.toWN();
        end
        
        function wn = multiplyMomentBased(this, wn2)
            % Multiply two WN distributions (approximate by means of true moments of product)
            %
            % Parameters:
            %   wn2 (WNDistribution)
            %       distribution to multiply with
            % Returns:
            %   wn (WNDistribution)
            %       product of this and wn2
            %
            % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
            % Recursive Bayesian Filtering in Circular State Spaces
            % arXiv preprint: Systems and Control (cs.SY), January 2015.
            arguments
                this (1,1) WNDistribution
                wn2 (1,1) WNDistribution
            end
            
            %integral(@(x) exp(1i*x).*wn.pdf(x).*wn.pdf(x), 0, 2*pi)/integral(@(x) wn.pdf(x).*wn.pdf(x), 0, 2*pi)
            erfz = @(z) Faddeeva_erf(z);
            
            sum1 = 0; %this sum is integral of exp(ix)*f1(x)*f2(x)
            sum2 = 0; %this sum is integral of f1(x)*f2(x)
            n = 10;
            for j=-n:n
                for k=-n:n
                    mu1 = this.mu + 2*pi*j;
                    mu2 = wn2.mu + 2*pi*k;
                    m = (mu1 * wn2.sigma^2 + mu2 * this.sigma^2)/(wn2.sigma^2+this.sigma^2);
                    s = sqrt((this.sigma^2*wn2.sigma^2)/(wn2.sigma^2+this.sigma^2));
                    sum1 = sum1 + exp(-1/2*(mu1-mu2)^2/(wn2.sigma^2+this.sigma^2)) * 1/2 * exp(1i*m - s^2/2) * (erfz((m + 1i*s^2)/(sqrt(2)*s)) - erfz((m - 2*pi + 1i*s^2)/(sqrt(2)*s)));
                    sum2 = sum2 + exp(-1/2*(mu1-mu2)^2/(wn2.sigma^2+this.sigma^2)) * 1/2 * (erf(m/(sqrt(2)*s)) - erf((m - 2 * pi)/(sqrt(2)*s)));
                end
            end
            moment = sum1/sum2;
            wn = WNDistribution.fromMoment(moment);
        end
        
        function wn = multiplyTruncated(this, wn2)
            % Multply two WN distributions (approximate using truncated series)
            %
            % Parameters:
            %   wn2 (WNDistribution)
            %       distribution to multiply with
            % Returns:
            %   wn (WNDistribution)
            %       product of this and wn2
            %
            % see Traa, J. Multichannel Source Separation and Tracking with
            % Phase Differences by Random Sample Consensus University of
            % Illinois, 2013 (chapter 4)
            %
            % Traa, J. & Smaragdis, P. A Wrapped Kalman Filter for Azimuthal
            % Speaker Tracking IEEE Signal Processing Letters, 2013, 20,
            % 1257-1260
            %
            % warning: this is not symmetric
            arguments
                this (1,1) WNDistribution
                wn2 (1,1) WNDistribution
            end
            
            n = 5; %influences number of summands
            K = this.sigma^2/(this.sigma^2 + wn2.sigma^2); %equivalent to Kalman gain
            g = 0;
            for l=-n:n
                eta = normpdf(wn2.mu+2*pi*l, this.mu, wn2.sigma)/wnpdf(wn2.mu+2*pi*l, this.mu, wn2.sigma);
                g = g + (wn2.mu + 2*pi*l - this.mu)*eta;
            end
            mu_ = mod(this.mu + K*g,2*pi);
            sigma_ = sqrt((1-K)*this.sigma^2);
            wn = WNDistribution(mu_, sigma_);
        end
        
        function wn = convolve(this, wn2)
            % Calculate convolution of two WN distributions (exact)
            %
            % Parameters:
            %   wn2 (WNDistribution)
            %       distribution to convolve with
            % Returns:
            %   wn (WNDistribution)
            %       convolution of this and wn2
            arguments
                this (1,1) WNDistribution
                wn2 (1,1) WNDistribution
            end
            
            mu_ = mod(this.mu+wn2.mu,2*pi);
            sigma_ = sqrt(this.sigma^2+wn2.sigma^2);
            wn = WNDistribution(mu_, sigma_);
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
            s = mod(this.mu + this.sigma * randn(1,n),2*pi);
        end
        
        function [samples, weights] = sampleOptimalQuantization(this, N)
            % Computes optimal quantization of the
            % wrapped normal distribution.
            %
            % Parameters:
            %   N (scalar)
            %       number of samples
            % Returns:
            %   samples (1 x N)
            %       N samples on the circle
            %   weights (1 x N)
            %       weight for each sample
            %
            % Igor Gilitschenski, Gerhard Kurz, Uwe D. Hanebeck, Roland Siegwart,
            % Optimal Quantization of Circular Distributions
            % Proceedings of the 19th International Conference on Information Fusion (Fusion 2016), Heidelberg, Germany, July 2016.
            
            assert(isscalar(N));
            
            funWithDeriv = @(x) optimfun(x, this.sigma);          
            startval = ((2*pi/(2*N) + (0:N-1)*(2*pi)/N)-pi) * min(this.sigma,1);
            
            [x, ~, ~, ~] = ...
                fminunc(funWithDeriv, startval, ...
                optimoptions( 'fminunc', 'Display','none', 'GradObj', 'on', ...
                'MaxFunEvals', 100000, 'MaxIter', 1000, ...
                'TolFun', 1e-12, 'TolX', 1e-12, 'Algorithm', 'trust-region'));
            
            samples = sort(mod(x + this.mu,2*pi));
            weights = WNWeights(samples, this.mu, this.sigma);
            
            function [v, g] = optimfun(x,sigma)
                v = WNUtility(x, sigma);
                g = WNUtilityDeriv(x, sigma);
            end
            
            function R = WNUtility( x, sigma )
                % Utility function for computing optimal L-Quantizer
                % Utility function used in optimization procedure computing the optimal
                % L-Quantizer of the wrapped normal distribution with mean 0
                %
                % Parameters:
                %       x - Positions of the discrete points (must be between -pi and pi).
                %       sigma - Dispersion parameter of wrapped normal.
                
                assert(sigma > 0, 'Sigma must be a positive real number.');
                
                R = sum(P1(x,sigma,0)-2*P2(x, sigma, 0)+P3(x, sigma, 0));
                i = 1;
                while 1
                    T =   P1(x,sigma,i) - 2*P2(x, sigma, i)  + P3(x, sigma, i) ...
                        + P1(x,sigma,-i) - 2*P2(x, sigma, -i) + P3(x, sigma, -i);                    
                    R = R + sum(T);
                    if (sum(T) == 0)
                        break;
                    end
                    i = i + 1;
                end
            end
            
            function res = P1(x, sigma, k)
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.
                
                xBa = xBoundary(1:(end-1));
                xBb = xBoundary(2:end);
                
                P = sigma * (2*k*pi-xBb).* exp( -((xBb + 2*k*pi).^2)/(2*sigma^2) ) / (sqrt(2*pi)) ...
                    - sigma *(2*k*pi-xBa).* exp( -((xBa + 2*k*pi).^2)/(2*sigma^2) ) / (sqrt(2*pi));
                P = P + ( (2*k^2*pi^2+sigma^2/2)*erf( (xBb + 2*k*pi)/(sqrt(2)*sigma) ) ...
                    - (2*k^2*pi^2+sigma^2/2)*erf( (xBa + 2*k*pi)/(sqrt(2)*sigma) ) );
                
                res = P;
            end
            
            function res = P2(x, sigma, k)
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.
                
                xBa = xBoundary(1:(end-1));
                xBb = xBoundary(2:end);
                
                P = sigma *   exp( -((xBb + 2*k*pi).^2)/(2*sigma^2) ) / (sqrt(2*pi)) ...
                    - sigma * exp( -((xBa + 2*k*pi).^2)/(2*sigma^2) ) / (sqrt(2*pi));
                P = P + ( k*pi*erf( (xBb + 2*k*pi)/(sqrt(2)*sigma) ) ...
                    - k*pi*erf( (xBa + 2*k*pi)/(sqrt(2)*sigma) ) );
                
                res = -x(2:(end-1)).*P;
            end
            
            function res = P3(x, sigma, k)
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.
                
                P = erf( (xBoundary(2:end)+2*k*pi)/(sqrt(2)*sigma) ) ...
                    - erf( (xBoundary(1:(end-1))+2*k*pi)/(sqrt(2)*sigma) );
                
                res = (x(2:(end-1)).^2).*P/2;
            end
            
            function G = WNUtilityDeriv( x, sigma )
                % Derivative of utility function for computing optimal L-Quantizer
                % Utility function used in optimization procedure computing the optimal
                % L-Quantizer of the wrapped normal distribution with mean 0
                %
                % Parameters:
                %       x - Positions of the discrete points (must be between -pi and pi).
                %       sigma - Dispersion parameter of wrapped normal.
                
                assert(sigma > 0, 'sigma must be a positive real number.');
                
                G = zeros(1,numel(x));
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.
                
                dist = WNDistribution(0,sigma);
                WNDensity = @(x) dist.pdf(x);
                
                % Compute the actual utility.
                for i=1:(numel(x)-2)
                    intfun = @(a) (-2*a+2*x(i+1)).*WNDensity(a);
                    G(i) = integral(intfun, xBoundary(i),xBoundary(i+1));
                end
            end
            
            function w = WNWeights(x, mu, sigma )
                % Computes probability weights of a given quantizer
                x = x-mu;
                w = PMass(x, sigma, 0);
                
                i=1;
                while 1
                    T = PMass(x, sigma, i)+PMass(x, sigma, -i);
                    if (sum(T) == 0)
                        break;
                    end
                    w = w + T;
                    i=i+1;
                end
                
                % Normalization to avoid numerical inaccuracy.
                w=w/sum(w);
            end
            
            function res = PMass(x, sigma, k)
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.                
                
                P = erf( (xBoundary(2:end)+2*k*pi)/(sqrt(2)*sigma) ) ...
                    - erf( (xBoundary(1:(end-1))+2*k*pi)/(sqrt(2)*sigma) );
                
                res = P/2;
            end
        end
        
        function wn = shift(this, angle)
            % Shift distribution by the given angle
            %
            % Parameters:
            %   shiftAngles (scalar)
            %       angle to shift by
            % Returns:
            %   wn (WNDistribution)
            %       shifted distribution
            assert(isscalar (angle));
            wn = WNDistribution(this.mu+angle,this.sigma); % Mod is done in constructor
        end
        
        function g = toGaussian(this)
            % Convert to 1D Gaussian distribution
            % this is a simple conversion that just keeps the parameters
            % for large sigma, better conversions are possible
            %
            % Returns:
            %   g (GaussianDistribution)
            %       Gaussian with same parameters
            g = GaussianDistribution(this.mu, this.sigma^2);
        end
    end
    
    methods (Static)
        function wn = fromMoment(m)
            % Obtain a WN distribution from a given first trigonometric moment
            %
            % Parameters:
            %   m (scalar)
            %       first trigonemetric moment (complex number)
            % Returns:
            %   wn (WNDistribution)
            %       WN distribution obtained by moment matching
            assert(isscalar(m));
            
            mu_ = mod(atan2(imag(m),real(m)),2*pi);
            sigma_ = sqrt(-2*log(abs(m)));
            wn = WNDistribution(mu_, sigma_);
        end
        
        function wn = mleJensen(samples)
            % Obtain WN distribution from samples using MLE
            % Parameters:
            %   samples (1 x n vector)
            %       array of samples in [0,2pi)
            % Returns:
            %   wn (WNDistribution)
            %       WNDistribution that maximizes approximated likelihood
            %
            % An approximation based on Jensen's inequality is used to make
            % the mle tractable
            %
            % see
            % Roy, Anandarup and Parui, Swapan K. and Roy, Utpal, "SWGMM: A
            % Semi-wrapped Gaussian Mixture Model for Clustering of Circular-linear Data"
            % Pattern Analysis and Applications, 2014, Appendix A
            
            assert(size(samples,1)==1);
            assert(size(samples,2)>=2);
            
            % use moment-based solution as starting point
            wnMomentbased = WDDistribution(samples).toWN();
            mualpha = wnMomentbased.mu;
            sigmaalpha = wnMomentbased.sigma;
            
            % weighting function
            alpha = @(w,s) normpdf(s+2*w*pi, mualpha, sigmaalpha)/sum(normpdf(s+2*(-10:1:10)'*pi, mualpha, sigmaalpha));
            
            n = size(samples,2);
            mu_ = 0;
            sigmasq = 0;
            
            % calculate mu
            for i=1:n
                for w=-10:10
                    mu_ = mu_ + alpha(w,samples(i))*(samples(i)+2*w*pi);
                end
            end
            mu_=mu_/n;
            
            % calculate sigma^2
            for i=1:n
                for w=-10:10
                    sigmasq = sigmasq + alpha(w,samples(i))*(samples(i)+2*w*pi - mu_)^2;
                end
            end
            sigmasq=sigmasq/n;
            
            wn = WNDistribution(mu_,sqrt(sigmasq));
        end
        
        function wn = mleNumerical(samples)
            % Obtain WN distribution from samples using MLE (based on numerical optimization)
            % Parameters:
            %   samples (1 x n vector)
            %       array of samples in [0,2pi)
            % Returns:
            %   wn (WNDistribution)
            %       WNDistribution that maximizes likelihood
            %
            
            % log likelihood function, x(1)=mu, x(2)=sigma
            function y = f(x)
                if x(2)>0 %ensure that sigma of resulting distribution is positive
                    w = WNDistribution(x(1),x(2));
                    y = -w.logLikelihood(samples);
                else
                    y = Inf;
                end
            end
            
            %start value
            wd = WDDistribution(samples);
            wn = wd.toWN();
            startValue = [wn.mu, wn.sigma];
            
            %perform optimization
            result = fminsearch(@f,startValue,optimset('display','none'));
            wn = WNDistribution(result(1), result(2));
        end
    end
end

