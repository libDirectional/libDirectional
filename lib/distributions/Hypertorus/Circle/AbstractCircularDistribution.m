classdef (Abstract) AbstractCircularDistribution < AbstractHypertoroidalDistribution
    % AbstractCircularDistribution abstract base class for distributions on
    % the circle S1, or equivalently SO(2). The circle is parameterized as
    % [0,2pi).
    
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
                val=1-this.integralNumerical(xa,startingPoint);
            else
                val=this.integralNumerical(startingPoint,xa);
            end
        end
        
        function s = sampleMetropolisHastings(this, n, proposal, startPoint, burnIn, skipping)
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
            arguments
                this (1,1) AbstractDistribution
                n (1,1) {mustBePositive,mustBeInteger}
                proposal (1,1) function_handle = @()NaN % Default proposals are set in inherting classes
                startPoint (:,1) double = this.meanDirection() % Default starting points are set in inherting classes
                burnIn (1,1) double = 10
                skipping (1,1) double = 5
            end
            if nargin(proposal)==0
                wn = WNDistribution.fromMoment(this.trigonometricMoment(1));
                wn.mu = 0;
                proposal = @(x) mod(x + wn.sample(1), 2*pi);
            end
            s = sampleMetropolisHastings@AbstractDistribution(this, n, proposal, startPoint, burnIn, skipping);
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
            p2 = stem3(cos(theta), sin(theta), ftheta, 'Marker', 'none');
            if h==1
                hold off
            end
            p = [p1 p2];
        end
         
        function h = plotCdf(this, startingPoint, resolution, varargin)
            % Create an appropriate plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to plot/surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            arguments
                this (1,1) AbstractCircularDistribution
                startingPoint (1,1) double = 0
                resolution (1,1) = 128
            end
            arguments (Repeating)
                varargin
            end
            theta = linspace(0,2*pi,resolution);
            ftheta = this.cdf(theta, startingPoint);
            h = plot(theta, ftheta,varargin{:});    
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
        
        function wd = toDirac5SuperPosition(this, lambda)
            % Convert to wrapped Dirac mixture using superposition method
            %
            % Parameters:
            %   lambda (scalar or 1 x n column vector)
            %       number of superpositions or vector containing lambda's
            % Returns:
            %   wd (WDDistribution)
            %       Dirac approximation
            %
            % Gerhard Kurz, Igor Gilitschenski, Roland Y. Siegwart, Uwe D. Hanebeck,
            % Methods for Deterministic Approximation of Circular Densities
            % Journal of Advances in Information Fusion, 11(2):138-156, December 2016.
            d = [];
            w = [];
            wCenter = 0;
            m1 = abs(this.trigonometricMoment(1));
            m2 = abs(this.trigonometricMoment(2));
            w5min = (4*m1^2 - 4*m1 - m2 + 1)/(4*m1 - m2 - 3);
            w5max = (2*m1^2-m2-1)/(4*m1 - m2 - 3);      
            if isscalar(lambda)
                % lambda defines the number of approximations to combine
                % use heuristic choice for lambda's
                q = lambda;
                lambdamax = 2*sqrt(2)-2;
                lambdamin = 2*q/(q-1)*w5min/(w5min-w5max) - lambdamax*(q+1)/(q-1);
                lambdamin = max(lambdamin, 0);
                lambda = lambdamin + (1:q)/(q) * (lambdamax-lambdamin);
            else
                % lambda is a column vector of the chosen lambda's
                assert(size(lambda,1) == 1);
                assert(all(lambda>=0));
                assert(all(lambda<=1));
            end
            for i=1:length(lambda)
                % generate Diracs for each lambda
                w5 = w5min + lambda(i) * (w5max-w5min);
                w1 = (1 - w5)/4;
                c1 = 2/(1-w5)*(m1-w5);
                c2 = 1/(1-w5)*(m2-w5)+1;
                x2 = (2*c1 + sqrt(4*c1^2 - 8* (c1^2-c2)))/4;
                x1 = c1-x2;
                phi1 = acos(x1);
                phi2 = acos(x2);
                diracs = [-phi1 -phi2 phi1 phi2 0];   
                weights = [w1 w1 w1 w1 w5];
                d = [d diracs(1:4)];
                w = [w weights(1:4)];    
                wCenter = wCenter + w5;
            end
            wAll = [w, wCenter];
            wd = WDDistribution([d 0]+this.circularMean, wAll/sum(wAll));
        end
        
        function wd = toDiracBT(this, n, fixCircularMean, fixFirstTrigonometrictMoment)
            % Convert to wrapped Dirac mixture using binary tree.
            %
            % Parameters:
            %   n (scalar)
            %       number of diracs
            %   fixCircularMean (boolean)
            %       corrects circular mean after generating the samples
            %   fixFirstTrigonometrictMoment (boolean)
            %       corrects the first trigonometric moment after
            %       generating the samples (requires a small numerical
            %       optimization)
            % Returns:
            %   wd (WDDistribution)
            %       Dirac approximation with n components
            %
            % Gerhard Kurz, Igor Gilitschenski, Roland Y. Siegwart, Uwe D. Hanebeck,
            % Methods for Deterministic Approximation of Circular Densities
            % Journal of Advances in Information Fusion, 11(2):138-156, December 2016.
            if nargin < 3   
                fixCircularMean = true;
            end
            if nargin < 4
                fixFirstTrigonometrictMoment = false;
            end
            if fixFirstTrigonometrictMoment && ~ fixCircularMean
                %it is not possible to fix the first trigonometric moment
                %without fixing the circular mean first
                fixCircularMean = true; 
            end
            
            function diracPos = approximateBT(left, right, n)
                %todo shift to mu=0?
                
                middle = (left+right)/2;
                % Base case: only one dirac remaining
                if n==1 
                    % Possibility 1: use middle for faster calculation
                    diracPos = middle;
                    % Possibility 2: use center of mass for higher accuracy
                    % diracPos = integral(@(x) x.*this.pdf(x), left, right)/this.cdf(right, left);% todo: bug if cdf is zero!
                    return
                end

                % Calculate integrals and distribute diracs proportionally
                leftInt = this.cdf(middle, left);
                rightInt = this.cdf(right, middle);
                nLeft  = floor(n*leftInt /(leftInt+rightInt));
                nRight = floor(n*rightInt/(leftInt+rightInt));

                % Check number of Dirac components
                if nLeft + nRight == n-1
                    if n*leftInt /(leftInt+rightInt) - nLeft > n*rightInt/(leftInt+rightInt) - nRight
                        nLeft = nLeft + 1;    
                    else
                        nRight = nRight + 1;
                    end
                elseif nLeft + nRight == n-2
                    % This should probably never happen except perhaps through
                    % numerical inaccuracies.
                    nLeft = nLeft + 1;
                    nRight = nRight + 1;
                    warning('nLeft + nRight == n-2');
                end

                % Recursive calls for left and right half
                diracPos = [];
                if nLeft >=1
                    diracPos = [diracPos approximateBT(left, middle, nLeft)];
                end
                if nRight>=1
                    diracPos = [diracPos approximateBT(middle, right, nRight)];
                end
            end
            
            d = approximateBT(0, 2*pi, n);
            wd = WDDistribution(d);
            
            % Moment correction
            if fixCircularMean
                % Correct circular mean
                meanError = this.circularMean() - wd.circularMean();
                d = d + meanError;
                wd = WDDistribution(d);
            end
            
            function wdScaled = scaledWd(c)
                % Scales the wd distribution around its mean by a factor of c
                m = wd.circularMean();
                d = (wd.d - m); % remove mean
                d = mod(d+pi, 2*pi)-pi; %move between -pi and pi
                d = d*c; % scale, note that scaling can alter the circular mean
                wdTemp = WDDistribution(d);
                wdScaled = WDDistribution(d - wdTemp.circularMean() + m); %add mean back in
            end
            
            function y = goalFunMomentError(c)
                wdTmp = scaledWd(c);
                y = abs(wdTmp.trigonometricMoment(1) - this.trigonometricMoment(1));
            end
            
            if fixFirstTrigonometrictMoment
                % Correct first trigonometric moment
                % Requires one-dimensional nonlinear optimization
                c = fminunc(@goalFunMomentError, 1, optimoptions(@fminunc, 'display', 'off', 'TolX', 1E-12, 'TolFun', 1E-12, 'algorithm','quasi-newton'));
                wd = scaledWd(c);
            end
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
        
        function result = l2distanceCdfNumerical(this, other, startingPoint)
            arguments
                this AbstractCircularDistribution
                other AbstractHypertoroidalDistribution
                startingPoint (1,1) double = 0
            end
            % Numerically calculates the L2 distance of the cdfs
            %
            % Parameters:
            %   other (AbstractHypertoroidalDistribution)
            %       distribution to compare with
            % Returns:
            %   result (scalar)
            %       test result
            startingPoint=mod(startingPoint,2*pi);
            assert(other.dim==1,'cdfs are only implemented for circular distributions.');
            if ~isa(other,'AbstractCircularDistribution')
                other = other.toCircular;
            elseif isa(other,'WDDistribution')
                % Swap distributions
                othertmp=this;
                this=other;
                other=othertmp;
            end
            if isa(this,'WDDistribution') % Splitting the interal up is usually faster for WDDistributions
                [d1ShiftedSorted,sorting]=sort(mod(this.d-startingPoint,2*pi)+startingPoint);
                wSorted=this.w(sorting);

                d1ShiftedSortedAppended=[startingPoint,d1ShiftedSorted,startingPoint+2*pi];
                wSortedAppended=[0,wSorted];
                wCum=cumsum(wSortedAppended);
                if isa(other,'WDDistribution') % In this special case, we can calculate the result very efficiently
                    d2Shifted=mod(other.d-startingPoint,2*pi)+startingPoint;
                    allBorders=[startingPoint,union(d1ShiftedSorted,d2Shifted,'sorted'),startingPoint+2*pi];
                    % Integrate by multiplying lenght of interval with
                    % squared diff of cdfs
                    result=sum((allBorders(2:end)-allBorders(1:end-1)).*...
                        (this.cdf(allBorders(1:end-1),startingPoint)-other.cdf(allBorders(1:end-1),startingPoint)).^2);
                elseif isa(other,'FourierDistribution')&&(strcmp(other.transformation,'sqrt')||strcmp(other.transformation,'identity')) % Additional speedup for Fourier
                    if strcmp(other.transformation,'sqrt') %transform to identity
                        fd=other.transformViaCoefficients('square',4*length(other.a)-3);
                    else
                        fd=other;
                    end
                    % For the cdf
                    c=fd.c;
                    c0=c((length(c)+1)/2);
                    cnew=fd.c./(1i*(-length(fd.b):length(fd.b))); 
                    cnew((length(cnew)+1)/2)=1/(2*pi);
                    fdInt=FourierDistribution.fromComplex(cnew,'identity');
                    
                    % For the integral of the cdf
                    cInt=fdInt.c;
                    c0int=real(cInt((length(cInt)+1)/2));
                    cintInt=fdInt.c./(1i*(-length(fdInt.b):length(fdInt.b))); % For all entries except c0
                    cintInt((length(cintInt)+1)/2)=1/(2*pi);
                    fdIntInt=FourierDistribution.fromComplex(cintInt,'identity'); 
                    fIntInt=@(r,l)fdIntInt.value(r)-fdIntInt.value(l)+(c0int)*(r-l);

                    % For the integral of the square of the cdf
                    cIntSquare=conv(cInt,cInt);
                    c0intSquare=cIntSquare((length(cIntSquare)+1)/2);
                    cintSquareInt=cIntSquare./(1i*(-2*length(fdInt.b):2*length(fdInt.b))); % For all entries except c0
                    cintSquareInt((length(cintSquareInt)+1)/2)=1/(2*pi);
                    fdIntSquareInt=FourierDistribution.fromComplex(cintSquareInt,'identity'); 
                    fIntSquareInt=@(r,l)fdIntSquareInt.value(r)-fdIntSquareInt.value(l)+(c0intSquare)*(r-l);
                
                    d1Ssa=d1ShiftedSortedAppended;
                    sp=startingPoint;
                    result=... % Calculate vectorized over the weights and boundaries
                        ...% Parts dependent on y
                        sum((d1Ssa(2:end)-d1Ssa(1:end-1)).*(wCum.^2+2*fdInt.value(sp)*wCum+fdInt.value(sp)^2 ...
                             +2*fdInt.value(sp)*c0*sp+2*c0*sp*wCum+c0^2*sp^2))...
                        ...% Parts dependent on y^2
                        +1/3*sum( (d1Ssa(2:end).^3-d1Ssa(1:end-1).^3)*c0^2)...
                        ...% Parts dependent on fy^2
                        +sum(fIntSquareInt(d1Ssa(2:end),d1Ssa(1:end-1)))...
                        ...% Part dependent on fy
                        +sum(fIntInt(d1Ssa(2:end),d1Ssa(1:end-1)).*(-2*fdInt.value(sp)-2*wCum-2*c0*sp))...
                        ...% Part dependent on y
                        +sum(0.5*(d1Ssa(2:end).^2-d1Ssa(1:end-1).^2).*(-2*c0*wCum-2*sp*c0^2-2*fdInt.value(sp).*c0));
                    for i=2:length(d1Ssa) % Iterate over all intervals for calculation of part depending on y*fy
                        l=d1Ssa(i-1);
                        r=d1Ssa(i);
                        % Calculate integral for the part dependeing on y*fy
                        kRange=-(numel(cInt)-1)/2:(numel(cInt)-1)/2;
                        intValyfySummands=real(cInt.*(exp(l*kRange*1i).*(-1+l*kRange*1i)./kRange.^2-exp(r*kRange*1i).*(-1+r*kRange*1i)./kRange.^2));
                        intValyfySummands((numel(cInt)+1)/2)=0;
                        intValyfy=sum(intValyfySummands)+cInt((numel(cInt)+1)/2)*(r^2/2 - l^2/2); % Hanlde part depending on c0int separately

                        result=real(result)+2*c0*real(intValyfy);
                    end
                elseif isa(other,'FIGDistribution')
                    if distribution.enforcePdfNonnegative
                        desiredTransformation = 'sqrt';
                    else
                        desiredTransformation = 'identity';
                    end
                    % We take the same number of coefficients as is done
                    % when evaluating the pdf
                    fdOther = FourierDistribution.fromDistribution(other,desiredTransformation);
                    result = this.l2distanceCdfNumerical(fdOther);
                else
                    result=sum(arrayfun(@(i)...
                        integral(@(x)(other.cdf(x,startingPoint)-wCum(i-1)).^2,...
                            d1ShiftedSortedAppended(i-1),d1ShiftedSortedAppended(i)),...
                        2:length(d1ShiftedSortedAppended)));
                end
            else
                result=integral(@(x) (this.cdf(x,startingPoint)-other.cdf(x,startingPoint)).^2, 0, 2*pi);
            end
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

