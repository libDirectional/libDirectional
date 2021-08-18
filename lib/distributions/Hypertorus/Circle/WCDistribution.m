classdef WCDistribution < AbstractCircularDistribution
    % Wrapped Cauchy distribution
    %
    % see Sreenivasa Rao Jammalamadaka and A. SenGupta, "Topics in Circular 
    % Statistics", 2001, Sec. 2.2.7, page 45f.
    
    properties
        mu
        gamma
    end
    
    methods
        function this = WCDistribution(mu_, gamma_)
            % Constructor
            this.mu = mod(mu_,2*pi);
            assert(gamma_>0);
            this.gamma = gamma_;
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
            p = 1/2/pi * sinh(this.gamma)./(cosh(this.gamma) - cos(xa-this.mu));
        end
        
        function p = cdf(this, xa, startingPoint)
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
            
            % see Wolfram Alpha
            % https://www.wolframalpha.com/input/?i=Integral%281%2F2%2Fpi*sinh%28g%29%2F%28cosh%28g%29-cos%28x-mu%29%29%2C+x%29
            
            assert(size(xa,1)==1);
            if nargin<=2 
                startingPoint = 0; 
            end   
            
            f = @(x) atan(coth(this.gamma/2) * tan((x-this.mu)/2))/pi;
            p = mod(f(xa) - f(startingPoint),1) ;
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
            m = exp(1i*n*this.mu - abs(n)*this.gamma);
        end
        
        function wc = convolve(this, wc2)
            % Calculate convolution of two WC distributions (exact)
            %
            % Parameters:
            %   wc2 (WCDistribution)
            %       distribution to convolve with
            % Returns:
            %   wc (WCDistribution)
            %       convolution of this and wc2
            %
            % see https://en.wikipedia.org/wiki/List_of_convolutions_of_probability_distributions#Continuous_distributions
            % Dwass, M. On the Convolution of Cauchy Distributions
            % The American Mathematical Monthly, Mathematical Association of America, 1985, 92, 55-57
            assert(isa(wc2, 'WCDistribution'));
            
            wc = WCDistribution(this.mu + wc2.mu, this.gamma + wc2.gamma);
        end       
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            result = log(2*pi * (1-exp(-2*this.gamma)));
        end
        
        function wc = shift(this,angle)
            % Shift distribution by the given angle
            %
            % Parameters:
            %   shiftAngles (scalar) 
            %       angle to shift by
            % Returns:
            %   hd (WCDistribution)
            %       shifted distribution
            assert(isscalar (angle));
            wc = WCDistribution(this.mu+angle,this.gamma); % Mod is done in constructor
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

            % based on inverting the cdf of a Cauchy distribution analytically
            y = rand(1,n); % uniform samples
            s = mod(this.mu + this.gamma * tan( pi * y-1/2 ), 2*pi); 
        end
    end
    
    methods (Static)
        function wc = fromMoment(m)
            % Obtain a WC distribution from a given first trigonometric moment
            %
            % Parameters:
            %   m (scalar)
            %       first trigonemetric moment (complex number)
            % Returns:
            %   wc (WCDistribution)
            %       WC distribution obtained by moment matching
            mu_ = mod(atan2(imag(m),real(m)),2*pi);
            gamma_ = - log(abs(m));
            wc = WCDistribution(mu_, gamma_);
        end
    end
end

