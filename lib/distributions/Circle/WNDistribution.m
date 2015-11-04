classdef WNDistribution < AbstractCircularDistribution
    % wrapped normal distribution
    %
    % see Sreenivasa Rao Jammalamadaka and A. SenGupta, "Topics in Circular 
    % Statistics", 2001, Sec. 2.2.6, page 44f.    
    
    properties
        mu
        sigma
    end
    
    methods
        function this = WNDistribution(mu_, sigma_)
            % Constructor
            assert(isscalar(mu_));
            this.mu = mod(mu_,2*pi);
            assert(sigma_>0);
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
            assert(size(xa,1)==1);
            p = wnpdf(xa, this.mu, this.sigma);
        end
        
        function val=cdf(this,xa,startingPoint,n)
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
            if nargin<=2 
                startingPoint = 0; 
            end
            if nargin<=3 %Use n for predefined number of runs to use vectorization
                n=10; 
            end
            startingPoint=mod(startingPoint,2*pi);
            xa=mod(xa,2*pi);

            ncdf=@(from,to)1/2*(erf((this.mu-from)/(sqrt(2)*this.sigma))-erf((this.mu-to)/(sqrt(2)*this.sigma)));
            val=ncdf(startingPoint,xa);
            for i=1:n
                val=val+ncdf(startingPoint+2*pi*i,xa+2*pi*i)+ncdf(startingPoint-2*pi*i,xa-2*pi*i);
            end
            val=(xa<startingPoint)+val; %1+int(a,b) if a>b, otherweise int(a,b)
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
            
            assert(isa(wn2, 'AbstractCircularDistribution'));
            
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
            
            assert(isa(wn2, 'WNDistribution'));
            
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
            
            assert(isa(wn2, 'WNDistribution'));
            
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
            
            assert(isa(wn2, 'WNDistribution'));
            
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
        
        function wn=shift(this,angle)
            wn=WNDistribution(this.mu+angle,this.sigma); % Mod is done in constructor
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

