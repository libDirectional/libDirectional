classdef WDDistribution < AbstractCircularDistribution & HypertoroidalWDDistribution
    % Wrapped Dirac distribution (i.e., discrete distribution)
    % with Dirac positions d and weights w
    %
    % see Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Recursive Nonlinear Filtering for Angular Data Based on Circular Distributions
    % Proceedings of the 2013 American Control Conference (ACC 2013), Washington D. C., USA, June 2013
    
    methods
        function this = WDDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (1 x L)
            %       Dirac locations in [0,2pi)
            %   w_ (1 x L)
            %       weights for each Dirac
            arguments
                d_ (1,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            this@HypertoroidalWDDistribution(d_, w_);
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
            % shift d and xa so that startingPoint is beginning of the 2pi
            % interval
            arguments
                this (1,1) WDDistribution
                xa (1,:) double
                startingPoint (1,1) double = 0
            end
            dshifted = this.d+(this.d<startingPoint)*2*pi;
            xamod = mod(xa,2*pi);
            xamodshift = xamod+(xamod<startingPoint)*2*pi; 
            % If we wish to evaluate the cdf at many xa at once, there are
            % multiple potential ways to influence the performance. The
            % function could be converted into a piecewise constant 
            % function (with borders that can be freely chosen) and we
            % could then apply binary search to find the suitable region.
            % However, without an efficient binary search directly
            % available in Matlab, we have not tried this. Naive
            % futher vectorization that results in replication of data, e.g.,
            % val=(((ones(length(xa),1)*dshifted)<=(xamodshift'*ones(1,length(this.d))))*this.w')';
            % has been tested but is significantly slower than the approach
            % implemented below.
            val = arrayfun(@(xi)sum(this.w(dshifted<=xi)),xamodshift);
        end
                
        function m = trigonometricMoment(this,n)
            arguments
                this (1,1) WDDistribution
                n (1,1) double {mustBeInteger}
            end
            % Calculate n-th trigonometric moment analytically
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th trigonometric moment (complex number)
            m = this.w*exp(1i*n*this.d).';
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
            arguments
                this (1,1) WDDistribution
                other (1,1) AbstractCircularDistribution
            end
            [dsorted, order] = sort(this.d);
            wsorted = this.w(order);
            z = other.cdf(dsorted);
            dplus = max(cumsum(wsorted) - z);
            dminus = max(z - [0, cumsum(wsorted(1:end-1))]);
            v = dplus + dminus;
        end
                        
        
        function pwc = toPWC(this, n)
            % Approximate the WD distribution using a piecewise constant
            % (PWC) distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of intervals to use
            % Returns:
            %   pwc (PWCDistribution)
            %       resulting approximation
            arguments
                this (1,1) WDDistribution
                n (1,1) double {mustBeInteger,mustBePositive}
            end
            wpwc = zeros(1,n);
            for i=1:n
                l = PWCDistribution.leftBorder(i,n);
                r = PWCDistribution.rightBorder(i,n);
                wpwc(i) = sum(this.w(this.d >= l & this.d < r));
            end
            pwc = PWCDistribution(wpwc);
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
            arguments
                this (1,1) WDDistribution
            end
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
                
        function cd = toContinuousVoronoi(this)
            % Converts the wd to a continouos function using an approach
            % whose idea comes from Voronoi-regions. Find the center
            % between each two particles and use these as a boundary to
            % distribute the probability mass on (the probability mass in
            % this region corresponds to the weight of the dirac). Subsume
            % weights diracs with identical position beforehand
            arguments
                this (1,1) WDDistribution
            end
            [diracPositions,order]=sort(this.d); %Sort particles
            weights=this.w(order);
            % Calculate difference to next sample and borders
            difftmp=diff([0,diracPositions,2*pi]);
            % Now calculate distance from each sample to its next (wrap
            % around at border of periodicity)
            differences=[difftmp(2:end-1),difftmp(1)+difftmp(end)];
            for i=length(diracPositions)-1:-1:1
                if differences(i)==0
                    weights(i)=weights(i)+weights(i+1);
                    diracPositions(i+1)=[];
                    weights(i+1)=[];
                end
            end
            % Subdivide into regions
            difftmp=diff([0,diracPositions,2*pi]);
            differences=[difftmp(2:end-1),difftmp(1)+difftmp(end)];
            
            boundaries=diracPositions+differences/2;
            if boundaries(end)<2*pi
                allBoundaries=[0,boundaries(1:end),2*pi];
                weights=[weights,weights(1)];
            else
                allBoundaries=[0,boundaries(end)-2*pi,boundaries(1:end-1),2*pi];
                weights=[weights(end),weights];
            end
            regionSizes=[allBoundaries(2)-(allBoundaries(end-1)-2*pi),diff(allBoundaries(2:end-1)),allBoundaries(2)-(allBoundaries(end-1)-2*pi)];
            densityEachRegion=weights./regionSizes;
            
            cd=CustomCircularDistribution(@(x)...
                arrayfun(@(xCurr)densityEachRegion((xCurr>=allBoundaries(1:end-1))&(xCurr<allBoundaries(2:end))),...
                mod(x,2*pi)));
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
            arguments
                this (1,1) WDDistribution
                l (1,1) double = 0
                r (1,1) double = 2*pi
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
            %   h (scalar)
            %       plot handle
            h = this.plot(varargin{:});
        end
        
        function h = plot2dcircular(this, varargin)
            % Create a circular 2D plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot command
            % Returns:
            %   h (scalar)
            %       plot handle
            h = line([cos(this.d); (1+this.w).*cos(this.d)], [sin(this.d); (1+this.w).*sin(this.d)]);
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
        
        function sampleCdf(~, ~)
            % Disable sampling algorithm relying on cdf
            error('PDF:UNDEFINED', 'not supported');
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
                for l=1:Lw
                    w_(1,j+L*(l-1)) = this.w(j) * wd2.w(l);
                    d_(1,j+L*(l-1)) = this.d(j) + wd2.d(l);
                end
            end
            wd = WDDistribution(d_,w_);
        end
    end   
    
    methods (Static)
        function wd = fromGridDistribution(dist)
            arguments
                dist (1,1) FIGDistribution
            end
            d = dist.grid';
            w = dist.gridValues'/sum(dist.gridValues);
            wd = WDDistribution(d,w);
        end
        
        function f = fromDistribution(distribution, noOfSamples)
            arguments
                distribution (1,1) AbstractCircularDistribution
                noOfSamples (1,1) {mustBePositive, mustBeInteger}
            end
            f = WDDistribution(distribution.sample(noOfSamples),ones(1,noOfSamples)/noOfSamples);
        end
    end
end

