classdef AbstractCircularDistribution < AbstractHypertoroidalDistribution
    % AbstractCircularDistribution abstract base class for distributions on
    % the circle S1, or equivalently SO(2). The circle is parameterized as
    % [0,2pi).
    
    methods (Abstract)
        % Evaluate pdf at positions stored in xa
        pdf(this, xa);
    end
    
    methods
        function this = AbstractCircularDistribution
            % Constructor
            %
            % This is a one-dimensional special case of the
            % AbstractHypertoroidalDistribution.
            this.dim = 1;
        end
        function val = cdfNumerical(this, xa, startingPoint)
            % Evaluate cumulative distribution function using numerical integration
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
            assert(size(xa,1) == 1),
            if nargin<=2
                startingPoint = 0;
            end
            if length(xa)>1
                val=arrayfun(@(xi)this.cdf(xi,startingPoint),xa);
                return
            end
            startingPoint=mod(startingPoint,2*pi);
            xa=mod(xa,2*pi);
            if xa<startingPoint
                val=1-this.integral(xa,startingPoint);
            else
                val=this.integral(startingPoint,xa);
            end
        end
        
        function val = cdf(this, xa, startingPoint)
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
            val = cdfNumerical(this, xa, startingPoint); 
        end
            
        function p = plot2d(this, varargin)
            % Create a 2D plot of the pdf (this function is for backwards compatibility)
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot command
            % Returns:
            %   p (scalar)
            %       plot handle
            p = this.plot(varargin{:});
        end
        
        function p = plot2dcircular(this, varargin)
            % Create a circular 2D plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot command
            % Returns:
            %   p (scalar)
            %       plot handle
            theta = linspace(0,2*pi,320);
            theta = [theta 0];
            ftheta = this.pdf(theta);
            scale = 1/max(ftheta);
            p = plot((1+scale*ftheta).*cos(theta),(1+scale*ftheta).*sin(theta),varargin{:});
        end
        
        function p = plot3d(this, varargin)
            % Create a 3D plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot command
            % Returns:
            %   p (1 x 2)
            %       plot handles for plot3 and stem3
            theta = [linspace(0,2*pi,128), 0];
            ftheta = this.pdf(theta);
            p1 = plot3(cos(theta), sin(theta), ftheta, varargin{:});
            if ~ishold
                hold on
                h=1;
            else
                h=0;
            end
            p2 = stem3(cos(theta), sin(theta), ftheta, 'Marker', 'none',varargin{:});
            if h==1
                hold off
            end
            p = [p1 p2];
        end
         
        function v = circularVariance(this)
            % Calculate the circular variance
            %
            % Returns:
            %   v (scalar)
            %       circular variance in [0, 1]
            a = this.trigonometricMoment(1);
            v = 1 - abs(a);
        end
        
        function wd = toDirac2(this)
            % Convert to wrapped Dirac mixture with 2 components
            %
            % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
            % Deterministic Approximation of Circular Densities with Symmetric 
            % Dirac Mixtures Based on Two Circular Moments
            % Proceedings of the 17th International Conference on 
            % Information Fusion (Fusion 2014), Salamanca, Spain, July 2014.
            %
            % Returns:
            %   wd (WDDistribution)
            %       Dirac approximation with two components
            m1 = this.trigonometricMoment(1);
            mu = mod(atan2(imag(m1),real(m1)),2*pi);
            alpha = acos(abs(m1)); 
            a1 = mod(mu-alpha,2*pi);
            a2 = mod(mu+alpha,2*pi);
            wd = WDDistribution([a1 a2]);
        end
        
        function wd = toDirac3(this)
            % Convert to wrapped Dirac mixture with 3 components
            %
            % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
            % Recursive Nonlinear Filtering for Angular Data Based on Circular Distributions
            % Proceedings of the 2013 American Control Conference (ACC 2013), 
            % Washington D. C., USA, June 2013.
            %
            % Returns:
            %   wd (WDDistribution)
            %       Dirac approximation with three components
            m1 = this.trigonometricMoment(1);
            mu = mod(atan2(imag(m1),real(m1)),2*pi);
            alpha = acos(3/2*abs(m1)-1/2); 
            a1 = mod(mu-alpha,2*pi);
            a2 = mu;
            a3 = mod(mu+alpha,2*pi);
            wd = WDDistribution([a1 a2 a3]);
        end
        
        function wd = toDirac5(this, lambda)
            % Convert to wrapped Dirac mixture with 5 components
            %
            % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
            % Deterministic Approximation of Circular Densities with Symmetric 
            % Dirac Mixtures Based on Two Circular Moments
            % Proceedings of the 17th International Conference on 
            % Information Fusion (Fusion 2014), Salamanca, Spain, July 2014.
            %
            % Parameters:
            %   lambda (scalar)
            %       weighting parameter in [0,1] for sample at circular mean
            % Returns:
            %   wd (WDDistribution)
            %       Dirac approximation with five components
            m1 = abs(this.trigonometricMoment(1));
            m2 = abs(this.trigonometricMoment(2));
            if nargin<2
                lambda = 0.5;
            else
                assert(lambda>=0 && lambda<=1);
            end
            w5min = (4*m1^2 - 4*m1 - m2 + 1)/(4*m1 - m2 - 3);
            w5max = (2*m1^2-m2-1)/(4*m1 - m2 - 3);
            w5 = w5min + lambda * (w5max-w5min);
            w1 = (1 - w5)/4;
            c1 = 2/(1-w5)*(m1-w5);
            c2 = 1/(1-w5)*(m2-w5)+1;
            x2 = (2*c1 + sqrt(4*c1^2 - 8* (c1^2-c2)))/4;
            x1 = c1-x2;
            phi1 = acos(x1);
            phi2 = acos(x2);
            diracs = [-phi1 phi1 -phi2 phi2 0];   
            weights = [w1 w1 w1 w1 w5];
            wd = WDDistribution(diracs+this.circularMean,weights);
        end
                
        function vm = toVM(this)
            % Convert to von Mises by trigonometric moment matching
            %
            % Returns:
            %   vm (VMDistribution)
            %       distribution with same first trigonometric moment
            vm = VMDistribution.fromMoment(this.trigonometricMoment(1));
        end
        
        function wc = toWC(this)
            % Convert to wrapped Cauchy by trigonometric moment matching
            %
            % Returns:
            %   wc (VMDistribution)
            %       distribution with same first trigonometric moment
            wc = WCDistribution.fromMoment(this.trigonometricMoment(1));
        end
        
        function wn = toWN(this)
            % Convert to wrapped normal by trigonometric moment matching
            %
            % Returns:
            %   wn (VMDistribution)
            %       distribution with same first trigonometric moment
            wn = WNDistribution.fromMoment(this.trigonometricMoment(1));
        end

        function result = integral(this, l, r)
            % Calculates the integral of the pdf from l to r analytically
            % if possible, fall back to numerical calculation by default
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
            result = this.integralNumerical(l,r);
        end
        
        function result = integralNumerical(this, l, r)
            % Calculates the integral of the pdf from l to r numerically
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
            result = integral(@(x) this.pdf(x), l, r);
        end
                
        function s = sampleCdf(this, n)
            % Sampling algorithm based on calculation of the inverse of the 
            % distribution function with numerical integration
            % (does not always work depending on shape of distribution)
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (1 x n)
            %       n samples on the circle
            s = zeros(1,n);
            r = rand(1,n);
            for i=1:n
                %F = @(x) integral(@(t) this.pdf(t), 0, x)-r(i);
                F = @(x) this.cdf(x)-r(i);
                s(i) = fsolve(F, pi/2, optimset('display', 'none'));
            end
        end
        
        function s = sampleMetropolisHastings(this, n)
            % Metropolis Hastings sampling algorithm
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (1 x n)
            %       n samples on the circle
            %
            % Hastings, W. K. 
            % Monte Carlo Sampling Methods Using Markov Chains and Their Applications 
            % Biometrika, 1970, 57, 97-109
            
            burnin = 10;
            skipping = 5;
            
            totalSamples = burnin+n*skipping;
            s = zeros(1, totalSamples);
            x = this.circularMean; %start with mean
            wn = WNDistribution.fromMoment(this.trigonometricMoment(1));
            wn.mu = 0;
            proposal = @(x) mod(x + wn.sample(1), 2*pi);
            i = 1;
            pdfx = this.pdf(x);
            while i <= totalSamples 
                xNew = proposal(x); %generate new sample
                pdfxNew = this.pdf(xNew);
                a = pdfxNew/pdfx;
                if a > 1
                    %keep sample
                    s(i) = xNew;
                    x = xNew;
                    pdfx = pdfxNew;
                    i = i+1;
                else
                    r = rand(1);
                    if a > r
                        %keep sample
                        s(i) = xNew;
                        x = xNew;
                        pdfx = pdfxNew;
                        i = i+1;
                    else
                        %reject sample
                    end
                end
            end
            s = s(burnin+1:skipping:end);
        end
    end

    methods(Static)
        function p = plotCircle(varargin)
            % Plots a circle, useful for 3D plot and circular 2D plot
            %
            % Returns:
            %   p (scalar)
            %       plot handle
            theta = [linspace(0,2*pi,320),0];
            p = plot(cos(theta),sin(theta),varargin{:});
        end
    end
end

