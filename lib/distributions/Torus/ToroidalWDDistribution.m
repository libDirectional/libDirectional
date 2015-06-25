classdef ToroidalWDDistribution < AbstractToroidalDistribution
    % Wrapped Dirac distribution with dirac positions d and
    % weights w.
    %
    % Gerhard Kurz, Igor Gilitschenski, Maxim Dolgov, Uwe D. Hanebeck,
    % Bivariate Angular Estimation Under Consideration of Dependencies Using Directional Statistics
    % Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014), Los Angeles, California, USA, December 2014.
    
    properties
        d
        w
    end
    
    methods
        function this = ToroidalWDDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (3 x L)
            %       Dirac locations in [0,2pi)^2
            %   w_ (1 x L)
            %       weights for each Dirac            
            assert(size(d_,1) ==2 );
            this.d = mod(d_, 2*pi);
            if (nargin<2)
                %all diracs have equal weights
                this.w = ones(1,size(this.d,2))/size(this.d,2);
            else
                assert(size(w_,1) == 1 );
                assert(size(d_,2) == size(w_,2));
                this.w = w_/sum(w_);
            end
        end
        
        function p = pdf(this, xa)
            % Placeholder, pdf does not exist for wrapped Dirac distributions
            p = 0; %ToroidalWDDistribution does not have a proper pdf
            warning('PDF:UNDEFINED', 'pdf is not defined')
        end
        
         function m = trigonometricMoment(this,n)
            % Calculate n-th trigonometric moment, i.e., 
            % E([e^(inx_1); e^(inx_2)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (2x1)
            %       n-th trigonometric moment (complex vector)  
            assert(isscalar(n));
            
            m(1,1) = sum(exp(1i*n*this.d(1,:)).*this.w);
            m(2,1) = sum(exp(1i*n*this.d(2,:)).*this.w);
         end
        
        function m = trigonometricMomentNumerical(this,n)
            % Disable numerical calculation of angular moments since it relies on the pdf
            error('PDF:UNDEFINED', 'not supported');
        end
        
        function m = angularProductMomentNumerical(this,n)
            % Disable numerical calculation of angular product moments since it relies on the pdf
            error('PDF:UNDEFINED', 'not supported');
        end
        
        function l = logLikelihood(this, samples)
            error('PDF:UNDEFINED', 'not supported');
        end           
        
        function pm = angularProductMoment(this, n)
            % Calculate E(e^inx * e^iny) as used by Jammalamadaka 1988
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th angular product moment (complex number)  
            assert(isscalar(n));
            
            m = this.circularMean();
            pm = sum(this.w.* exp(1i*n*(this.d(1,:)-m(1))).* exp(1i*n*(this.d(2,:)-m(2))));
        end
              
        function twd = applyFunction(this,f)
            % Apply a function f(x) to each Dirac component and obtain its new position
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi)^2 to [0,2*pi)^2
            % Returns:
            %   twd (ToroidalWDDistribution)
            %       distribution with new Dirac locations (and same
            %       weights as before)
            assert(isa(f,'function_handle'));
            
            d_ = zeros(size(this.d));
            for i=1:size(this.d,2)
                d_(:,i) = f(this.d(:,i));
            end
            twd = ToroidalWDDistribution(d_, this.w);
        end
        
        function twd = reweigh(this, f)
            % Uses a function f(x) to calculate the weight of each Dirac
            % component. The new weight is given by the product of the old 
            % weight and the weight obtained with f. Restores normalization
            % afterwards.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi)^2 to [0, infinity)
            % Returns:
            %   twd (ToroidalWDDistribution)
            %       distribution with new weights and same Dirac locations
            assert(isa(f,'function_handle'));
            
            wNew = zeros(1, length(this.w));
            for i=1:length(this.d)
                wNew(i) = f(this.d(:,i));
            end
            twd = ToroidalWDDistribution(this.d, wNew .* this.w);
        end
        
        function result = integralNumerical(this, l1, r1, l2, r2)
            error('PDF:UNDEFINED', 'not supported')
        end
        
        function covariance4DNumerical(this)
            error('PDF:UNDEFINED', 'not supported')
        end
        
        function result = integral(this, l1, r1, l2, r2)
            % Calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l1 (scalar)
            %       left bound of integral in first dimension, default 0
            %   r1 (scalar)
            %       right bound of integral in first dimension, default 2*pi          
            %   l2 (scalar)
            %       left bound of integral in second dimension, default 0
            %   r2 (scalar)
            %       right bound of integral in second dimension, default 2*pi          
            % Returns:
            %   result (scalar)
            %       value of the integral
            
            if nargin < 2;  l1 = 0;  end
            if nargin < 3;  r1 = 2*pi; end
            if nargin < 4;  l2 = 0; end
            if nargin < 5;  r2 = 2*pi; end
            
            % todo handle case where [l,r] spans more than 2pi
            assert(l1>=0 && l1<=2*pi);
            assert(r1>=0 && r1<=2*pi);
            assert(l2>=0 && l2<=2*pi);
            assert(r2>=0 && r2<=2*pi);
            
            if r1 < l1
               result = -this.integral(r1,l1, l2, r2); %swap l1 and r1
            elseif r2 < l2 
                result = -this.integral(l1,r1, r2, l2); %swap l2 and r2
            else
                %now we can guarantee l1<=r1 and l2<=r2
                result =  sum(this.w(this.d(1,:)>=l1 & this.d(1,:)<=r1 & this.d(2,:)>=l2 & this.d(2,:)<=r2 ));
            end
        end
        
        function h = plot(this, varargin)
            % Create a R^2->R plot 
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to stem command
            % Returns:
            %   p (scalar)
            %       plot handle
            h = stem3(this.d(1,:),this.d(2,:), this.w, varargin{:});
        end
        
        function p = plotCylinder(this, varargin)
            % Create a cylinder plot 
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to scatter3 command
            % Returns:
            %   p (scalar)
            %       plot handle
            offset = 0.5;
            scale = 0;
            p = scatter3(this.d(1,:), (offset+scale*this.w).*cos(this.d(2,:)), (offset+scale*this.w).*sin(this.d(2,:)), varargin{:});
        end
        
        function p = plotTorus(this, varargin)
            % Create a torus plot 
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to scatter3 command
            % Returns:
            %   p (scalar)
            %       plot handle
            scale = 0;% 1/max(max(f));
            a = 2; %larger radius
            b = 0.5; %smaller radius
            X = (a+(b+this.w*scale).*cos(this.d(1,:))).*cos(this.d(2,:));
            Y = (a+(b+this.w*scale).*cos(this.d(1,:))).*sin(this.d(2,:));
            Z = b.*sin(this.d(1,:));
            p = scatter3(X,Y,Z, varargin{:});
        end

        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (2 x n)
            %       one sample per column
            assert(isscalar(n));
            ids = discretesample(this.w,n);
            s = this.d(:,ids);
        end
               
        function s = sampleMetropolisHastings(this, n)
            % Disable sampling algorithm relying on pdf
            error('PDF:UNDEFINED', 'not supported');
        end
        
        function rhoc = circularCorrelationJammalamadaka(this)
            % Calculates Jammalamadaka's correlation coefficient
            % see Jammalamadaka, S. R. & Sarma, Y. 
            % A Correlation Coefficient for Angular Variables 
            % Statistical Theory and Data Analysis II, 1988, 349-364
            %
            % Returns:
            %   rhoc (scalar)
            %       Jammalamadaka's correlation coefficient     
            m = this.circularMean();
            x = sum( this.w .* sin(this.d(1,:)-m(1)) .* sin(this.d(2,:)-m(2)) );
            y = sqrt(sum(this.w .* sin(this.d(1,:)-m(1)).^2) * sum(this.w .* sin(this.d(2,:)-m(2)).^2 ));
            rhoc = x/y;
        end
        
        function twn = toToroidalWNjammalamadaka(this)
            % Gets the appropriate Torus WN by moment matching
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       TWN with same trigonometric moment and same E(e^inx * e^iny)
            %
            % see Jammalamadaka 2001, page 181
            % be aware that the formulas given there are for mu = [0;0] only!
            % note that the resulting C may not positive definite, so the conversion may fail 
            mu = this.circularMean();
            xbar = sum( this.w.* cos(this.d(1,:) - mu(1)));
            ybar = sum( this.w.* cos(this.d(2,:) - mu(2)));
            %zbar = sum( this.w.* cos( this.d(1,:) - mu(1) + this.d(2,:) - mu(2) ));
            %note that the imaginary part of zbar is ignored!
            zbar = real(sum( this.w.* exp(1i*(this.d(1,:) - mu(1) + this.d(2,:) - mu(2) ))));
            si1sqared = -2 * log(xbar);
            si2sqared = -2 * log(ybar);
            si12 = log(xbar*ybar/zbar);
            C = [si1sqared si12; si12 si2sqared];
            twn = ToroidalWNDistribution(mu,C);
        end
        
        function twn = toToroidalWN(this)
            % Gets the appropriate ToroidalWNDistribution by moment matching
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       TWN with same trigonometric moment and same correlation
            %       coefficient (Jammalamadaka's coefficient is used)
            %            
            % Note that the resulting C may not positive definite, so the conversion may fail 
            %
            % Gerhard Kurz, Igor Gilitschenski, Maxim Dolgov, Uwe D. Hanebeck,
            % Bivariate Angular Estimation Under Consideration of Dependencies Using Directional Statistics
            % Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014), Los Angeles, California, USA, December 2014.
            %
            % get si1squared and si2squared just as Jammalamadaka
            mu = this.circularMean();
            xbar = sum( this.w.* cos(this.d(1,:) - mu(1)));
            ybar = sum( this.w.* cos(this.d(2,:) - mu(2)));
            si1sqared = -2 * log(xbar);
            si2sqared = -2 * log(ybar);
            % obtain si12 by solving equation that matches correlation coefficient 
            a = sum( this.w.* sin(this.d(1,:) - mu(1)).*sin(this.d(2,:) - mu(2)));
            b = sum( this.w.* sin(this.d(1,:) - mu(1)).^2);
            c = sum( this.w.* sin(this.d(2,:) - mu(2)).^2);
            % a/sqrt(b*c) is the sample correlation coefficient
            z = sqrt(sinh(si1sqared)*sinh(si2sqared))* a/sqrt(b*c);
            % function for inverse of sinh
            arsinh = @(z) log(z + sqrt(z^2+1));
            si12 = arsinh(z);
            C = [si1sqared si12; si12 si2sqared];
            twn = ToroidalWNDistribution(mu,C);
        end
        
        function C = covariance4D(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), cos(x2), sin(x2)] 
            %
            % Returns:
            %   C (4 x 4)
            %       covariance matrix of [cos(x1), sin(x1), cos(x2), sin(x2)] 
            dbar = [cos(this.d(1,:))', sin(this.d(1,:))', cos(this.d(2,:))', sin(this.d(2,:))'];
            mu = this.w*dbar;
            n = length(this.d);
            C = (dbar-repmat(mu,n,1))'*diag(this.w)*(dbar-repmat(mu,n,1));
        end
        
        function wd = marginal(this,dimension)
            % Get marginal distribution in first or second dimension, i.e., 
            % f(x_1) or f(x_2), respectively
            %
            % Parameters:
            %   dimension (scalar)
            %       the marginal in which dimension to calculate (1 or 2),
            %       the other dimension is marginalized out
            % Returns:
            %   wd (WDDistribution)
            %       marginal distribution (marginals are WD-distributed)
            assert(dimension == 1 || dimension == 2);
            wd = WDDistribution(this.d(dimension,:), this.w);
        end
    end
    
end

