classdef ToroidalWDDistribution < AbstractToroidalDistribution & HypertoroidalWDDistribution
    % Wrapped Dirac distribution on the torus with dirac positions d and
    % weights w.
    %
    % Gerhard Kurz, Igor Gilitschenski, Maxim Dolgov, Uwe D. Hanebeck,
    % Bivariate Angular Estimation Under Consideration of Dependencies Using Directional Statistics
    % Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014), Los Angeles, California, USA, December 2014.
    
    methods
        function this = ToroidalWDDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (2 x L)
            %       Dirac locations in [0,2pi)^2
            %   w_ (1 x L)
            %       weights for each Dirac  
            arguments
                d_ (2,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            this@HypertoroidalWDDistribution(d_, w_);
        end
        
        function angularProductMomentNumerical(~,~)
            % Disable numerical calculation of angular product moments since it relies on the pdf
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
                
        function covariance4DNumerical(~)
            error('PDF:UNDEFINED', 'not supported')
        end
        
        function result = integral(this, l, r)
            % Calculates the integral of the pdf from l to r
            %
            % Parameters:
            %   l (2 x 1 column vector)
            %       left bound of integral in each dimension, default 0
            %   r (2 x 1 column vector)
            %       right bound of integral in each dimension, default 2*pi           
            % Returns:
            %   result (scalar)
            %       value of the integral
            if nargin < 2;  l = [0; 0]; end
            if nargin < 3;  r = [2*pi; 2*pi]; end
            assert(all(size(l) == [this.dim, 1]));
            assert(all(size(r) == [this.dim, 1]));
            
            % todo handle case where [l,r] spans more than 2pi
            assert(l(1)>=0 && l(1)<=2*pi);
            assert(r(1)>=0 && r(1)<=2*pi);
            assert(l(2)>=0 && l(2)<=2*pi);
            assert(r(2)>=0 && r(2)<=2*pi);
            
            if r(1) < l(1)
                result = -this.integral([r(1); l(2)], [l(1); r(2)]); %swap l1 and r1
            elseif r(2) < l(2) 
                result = -this.integral([l(1); r(2)], [r(1); l(2)]); %swap l2 and r2
            else
                %now we can guarantee l1<=r1 and l2<=r2
                result =  sum(this.w(this.d(1,:)>=l(1) & this.d(1,:)<r(1) & this.d(2,:)>=l(2) & this.d(2,:)<r(2) ));
            end
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
            scale = 0;% 1/max(f,[],[1,2]);
            a = 2; %larger radius
            b = 0.5; %smaller radius
            X = (a+(b+this.w*scale).*cos(this.d(1,:))).*cos(this.d(2,:));
            Y = (a+(b+this.w*scale).*cos(this.d(1,:))).*sin(this.d(2,:));
            Z = (b+this.w*scale).*sin(this.d(1,:));
            p = scatter3(X,Y,Z,this.w*10000, varargin{:});
            holdStatus = ishold;
            hold on
            this.plotTorusContours;
            if ~holdStatus
                hold off
            end
        end
        
        function p = plotTorusStemlike(this, varargin)
            % Create a torus plot 
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to scatter3 command
            % Returns:
            %   p (scalar)
            %       plot handle
            
            a = 2; %larger radius
            b = 0.5; %smaller radius
            
            scale = 0;
            X = (a+(b+this.w*scale).*cos(this.d(1,:))).*cos(this.d(2,:));
            Y = (a+(b+this.w*scale).*cos(this.d(1,:))).*sin(this.d(2,:));
            Z = b.*sin(this.d(1,:));
            
            scale = 1/max(this.w);
            X2 = (a+(b+this.w*scale).*cos(this.d(1,:))).*cos(this.d(2,:));
            Y2 = (a+(b+this.w*scale).*cos(this.d(1,:))).*sin(this.d(2,:));
            Z2 = (b+this.w*scale).*sin(this.d(1,:));
            
            for i=1:length(this.w)
                line([X(i); X2(i)], [Y(i); Y2(i)], [Z(i); Z2(i)], 'color', 'r');
            end
            p = scatter3(X2,Y2,Z2, varargin{:});
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
            % and page 357 of:
            % Jammalamadaka, S. R. & Sarma, Y. 
            % A Correlation Coefficient for Angular Variables 
            % Statistical Theory and Data Analysis II, 1988, 349-364
            %
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
            % Gets the appropriate ToroidalWNDistribution by moment
            % matching using Jammalamadaka's correlation coefficient
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
        
        function twn = toToroidalWNjupp(this)
            % Gets the appropriate ToroidalWNDistribution by moment
            % matching using Jupp's correlation coefficient
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       TWN with same trigonometric moment and same correlation
            %       coefficient (Jammalamadaka's coefficient is used)
            %            
            % Note that the resulting C may not positive definite, so the conversion may fail 
            %
            % Gerhard Kurz, Uwe D. Hanebeck,
            % Parameter Estimation for the Bivariate Wrapped Normal Distribution
            % Proceedings of the 54th IEEE Conference on Decision and Control (CDC 2015), Osaka, Japan, December 2015.
            %
            % use standard moment matching for everything except correlation
            mu = this.circularMean();
            xbar = sum( this.w.* cos(this.d(1,:) - mu(1)));
            ybar = sum( this.w.* cos(this.d(2,:) - mu(2)));
            c11 = -2 * log(xbar);
            c22 = -2 * log(ybar);

            rho = this.circularCorrelationJupp();

            % create preliminary twn to obtain certain entries of second moment
            C = [c11 0; 0 c22];
            twn = ToroidalWNDistribution([0;0],C);  
            C_ = twn.covariance4D();

            % obtain correlation
            eta = exp(-c11/2-c22/2);

            a = C_(2,2)*C_(4,4) + C_(1,1)*C_(3,3);
            b = - 2* C_(2,2)*C_(4,4);
            c =  C_(2,2)*C_(4,4) - C_(1,1)*C_(3,3) - rho^2 * C_(1,1)*C_(2,2)*C_(3,3)*C_(4,4) / eta^2 ;
            tau = (-b + sqrt(b^2-4*a*c))/(2*a);
            arccosh = @(x) log(x + sqrt(x^2-1));
            c12 = sign(rho)*arccosh(tau);
            
            %w1 = C_(1,1)*C_(3,3);
            %w2 = C_(2,2)*C_(4,4);
            %tau = (w2 + sqrt(w2^2 - (w2+w1)*(w2-w1-rho^2*w1*w2/eta^2)))/(w1+w2);
            %c12 = sign(rho)*arccosh(tau);

            % create final twn
            C = [c11 c12; c12 c22];
            twn = ToroidalWNDistribution(mu,C);
        end
        
        function twn = toToroidalWNcovariance(twd)
            % Gets the appropriate ToroidalWNDistribution by moment
            % matching using minimum of the Frobenius norm of the upper right 2x2
            % submatrix of the covariance
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       TWN with same trigonometric moment and same correlation
            %       coefficient (Jammalamadaka's coefficient is used)
            %            
            % Note that the resulting C may not positive definite, so the conversion may fail 
            %
            % Gerhard Kurz, Uwe D. Hanebeck,
            % Parameter Estimation for the Bivariate Wrapped Normal Distribution
            % Proceedings of the 54th IEEE Conference on Decision and Control (CDC 2015), Osaka, Japan, December 2015.
            %
            % use standard moment matching for everything except correlation
            mu = twd.circularMean();
            xbar = sum( twd.w.* cos(twd.d(1,:) - mu(1)));
            ybar = sum( twd.w.* cos(twd.d(2,:) - mu(2)));
            c11 = -2 * log(xbar);
            c22 = -2 * log(ybar);

            C_ = twd.covariance4D();
            eta = exp(-c11/2-c22/2);
            covss = C_(2,4);
            covcc = C_(1,3);

            c12candidates = ...
                [log(-(-covcc - covss - eta) / eta / 0.4e1 + sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) / 0.2e1 + sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.2e1 - 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 - (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1) + (-0.2e1 * (covcc - covss + eta) / eta - (-covcc - covss - eta) ^ 3 / eta ^ 3 / 0.4e1) * ((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) ^ (-0.1e1 / 0.2e1)) / 0.2e1)
                 log(-(-covcc - covss - eta) / eta / 0.4e1 + sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) / 0.2e1 - sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.2e1 - 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 - (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1) + (-0.2e1 * (covcc - covss + eta) / eta - (-covcc - covss - eta) ^ 3 / eta ^ 3 / 0.4e1) * ((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) ^ (-0.1e1 / 0.2e1)) / 0.2e1)
                 log(-(-covcc - covss - eta) / eta / 0.4e1 - sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) / 0.2e1 + sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.2e1 - 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 - (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1) - (-0.2e1 * (covcc - covss + eta) / eta - (-covcc - covss - eta) ^ 3 / eta ^ 3 / 0.4e1) * ((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) ^ (-0.1e1 / 0.2e1)) / 0.2e1)
                 log(-(-covcc - covss - eta) / eta / 0.4e1 - sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) / 0.2e1 - sqrt((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.2e1 - 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 - (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1) - (-0.2e1 * (covcc - covss + eta) / eta - (-covcc - covss - eta) ^ 3 / eta ^ 3 / 0.4e1) * ((-covcc - covss - eta) ^ 2 / eta ^ 2 / 0.4e1 + 0.1e1 / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (0.1e1 / 0.3e1) / 0.3e1 + (covcc ^ 2 + 0.2e1 * covcc * eta - covss ^ 2 - 0.3e1 * eta ^ 2) / eta * (-0.54e2 * covcc * covss * eta - 0.54e2 * eta ^ 2 * covss + 0.3e1 * sqrt(-0.3e1 * covcc ^ 6 - 0.18e2 * covcc ^ 5 * eta + 0.9e1 * covcc ^ 4 * covss ^ 2 - 0.9e1 * covcc ^ 4 * eta ^ 2 + 0.36e2 * covcc ^ 3 * covss ^ 2 * eta + 0.84e2 * covcc ^ 3 * eta ^ 3 - 0.9e1 * covcc ^ 2 * covss ^ 4 + 0.306e3 * covcc ^ 2 * covss ^ 2 * eta ^ 2 + 0.27e2 * covcc ^ 2 * eta ^ 4 - 0.18e2 * covcc * covss ^ 4 * eta + 0.540e3 * covcc * covss ^ 2 * eta ^ 3 - 0.162e3 * covcc * eta ^ 5 + 0.3e1 * covss ^ 6 + 0.27e2 * covss ^ 4 * eta ^ 2 + 0.405e3 * covss ^ 2 * eta ^ 4 + 0.81e2 * eta ^ 6)) ^ (-0.1e1 / 0.3e1)) ^ (-0.1e1 / 0.2e1)) / 0.2e1)];

            [~,best] = min(abs(imag(c12candidates))); % find solution with smallest imaginary part
            %todo check result
            c12 = real(c12candidates(best));
            C = [c11 c12; c12 c22];
            twn = ToroidalWNDistribution(mu,C);
        end
        
        function twn = toToroidalWNmixedMLE(this)
            % Gets the appropriate ToroidalWNDistribution by moment
            % matching for everything except the correlation and MLE for
            % the correlation
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       TWN with same trigonometric moment and same correlation
            %       coefficient (Jammalamadaka's coefficient is used)
            %            
            % This conversion should always succeed.
            %
            % Gerhard Kurz, Uwe D. Hanebeck,
            % Parameter Estimation for the Bivariate Wrapped Normal Distribution
            % Proceedings of the 54th IEEE Conference on Decision and Control (CDC 2015), Osaka, Japan, December 2015.
            %
            % use standard moment matching for everything except correlation
            mu = this.circularMean();
            xbar = sum( this.w.* cos(this.d(1,:) - mu(1)));
            ybar = sum( this.w.* cos(this.d(2,:) - mu(2)));
            si1sqared = -2 * log(xbar);
            si2sqared = -2 * log(ybar);
            C = [si1sqared 0; 0 si2sqared];
            twn = ToroidalWNDistribution(mu,C);

            function loglikelihood = calculateLogLikelihood(c12)
                twn.C(1,2) = c12;
                twn.C(2,1) = c12;
                if det(twn.C) > 0
                    loglikelihood = twn.logLikelihood(this.d, this.w);
                else
                    loglikelihood = -Inf;
                end
            end

            %c12min = -sqrt(twn.C(1,1))*sqrt(twn.C(2,2));
            %c12max = sqrt(twn.C(1,1))*sqrt(twn.C(2,2));

            %c12opt = fmincon(@(x) -calculateLogLikelihood(x), 0, [], [],[],[],c12min,c12max,[],optimset('display','none'));
            c12opt = fminsearch(@(x) -calculateLogLikelihood(x), 0, optimset('display','none'));

            twn.C(1,2) = c12opt;
            twn.C(2,1) = c12opt;
        end
        
        function twn = toToroidalWNunwrappingEM(this)
            % Gets the appropriate ToroidalWNDistribution by unwrapping
            % with the EM algorithm
            %
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       TWN with same trigonometric moment and same correlation
            %       coefficient (Jammalamadaka's coefficient is used)
            %            
            % This conversion should always succeed.
            %
            % "Time Series Analysis of Circular Data", N. I. Fisher and A.
            % J. Lee, Journal of the Royal Statistical Society. Series B (Methodological), 
            % Vol. 56, No. 2(1994), pp. 327-339
            %
            kmax = 1;

            epsilon = 1E-5;  % stopping condition
            maxIter = 100;% maximum number of iterations

            %start value
            mu=[0,0];
            C=eye(2,2);

            for iteration=1:maxIter
                samplesX=zeros(length(this.d),2*kmax+1,2*kmax+1);
                samplesY=zeros(length(this.d),2*kmax+1,2*kmax+1);
                p=zeros(length(this.d),2*kmax+1,2*kmax+1);
                
                % E-Step
                for i=1:length(this.d)
                    for k1=-kmax:kmax
                        for k2=-kmax:kmax
                            % probability that sample i was wrapped k1 times and k2
                            samplesX(i,k1+kmax+1,k2+kmax+1) = this.d(1,i)+ 2*pi*k1;
                            samplesY(i,k1+kmax+1,k2+kmax+1) = this.d(2,i)+ 2*pi*k2;
                            p(i,k1+kmax+1, k2+kmax+1) = this.w(i)*mvnpdffast([samplesX(i,k1+kmax+1,k2+kmax+1),samplesY(i,k1+kmax+1,k2+kmax+1)], mu, C);
                        end
                    end
                    % normalize per sample
                    p(i,:,:) =  p(i,:,:)/sum(p(i,:,:),[2,3]);
                end

                % M-Step
                p = p(:)/sum(p(:)); % normalize all
                muX = sum(p.*samplesX(:));
                muY = sum(p.*samplesY(:));
                muNew = [muX, muY];
                Cnew = [samplesX(:)-muX, samplesY(:)-muY]'*diag(p)*[samplesX(:)-muX, samplesY(:)-muY];
                
                % ensure symmetry
                Cnew = 1/2*(Cnew+Cnew');
                
                % check stopping condition
                if norm(muNew - mu) < epsilon && norm(Cnew(:) - C(:)) < epsilon
                    mu = muNew;
                    C = Cnew;    
                    break;
                end
                
                mu = muNew;
                C = Cnew;
            end

            twn = ToroidalWNDistribution(mu',C);
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
    end
    
end

