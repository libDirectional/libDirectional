classdef AbstractToroidalDistribution
    % Abstract base class for distributions on the torus (S1 x S1)
    % The torus is parameterized as [0,2*pi)^2
    
    properties
        
    end
    
    methods (Abstract)
        % Evaluate pdf at positions stored in xa
        pdf(this, xa);
    end
    
    methods
        function p = plot(this, varargin)
            % Create a surface plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            step = 2*pi/100;
            [alpha,beta] = meshgrid(0:step:2*pi,0:step:2*pi);
            f = this.pdf([alpha(:)'; beta(:)']);
            f = reshape(f,size(alpha,1), size(alpha,2));
            p = surf(alpha, beta, f, varargin{:});
        end
        
        function p = plotCylinder(this, varargin)
            % Create a cylinder plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            step = 2*pi/100;
            [alpha,beta] = meshgrid(0:step:2*pi,0:step:2*pi+step);
            f = this.pdf([alpha(:)'; beta(:)']);
            f = reshape(f,size(alpha,1), size(alpha,2));
            scale = 1/max(max(f));
            offset = 0.5;
            p = surf(alpha, (offset+scale*f).*cos(beta), (offset+scale*f).*sin(beta), f, varargin{:});
        end
        
        function p = plotTorus(this, varargin)
            % Create a torus plot of the pdf
            %
            % Parameters:
            %   varargin
            %       parameters to be passed to surf command
            % Returns:
            %   p (scalar)
            %       plot handle
            step = 2*pi/100;
            [alpha,beta] = meshgrid(0:step:2*pi,0:step:2*pi+step);
            f = this.pdf([alpha(:)'; beta(:)']);
            f = reshape(f,size(alpha,1), size(alpha,2));
            scale = 0;% 1/max(max(f));
            a = 2; %larger radius
            b = 0.5; %smaller radius
            X = (a+(b+f*scale).*cos(alpha)).*cos(beta);
            Y = (a+(b+f*scale).*cos(alpha)).*sin(beta);
            Z = b.*sin(alpha);
            p = surf(X,Y,Z,f, varargin{:});
        end
        
        function m = trigonometricMoment(this, n)
            % Calculate n-th trigonometric moment, i.e., 
            % E([e^(inx_1); e^(inx_2)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (2x1)
            %       n-th trigonometric moment (complex vector)             
            m = trigonometricMomentNumerical(this, n);
        end
        
        function m = trigonometricMomentNumerical(this, n)
            % Calculate n-th trigonometric moment numerically, i.e., 
            % E([e^(inx_1); e^(inx_2)])
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (2x1)
            %       n-th trigonometric moment (complex vector)             
            f1 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*exp(1i*n*x(:)'), size(x,1),size(x,2));
            f2 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*exp(1i*n*y(:)'), size(x,1),size(x,2));
            m(1,1) = integral2(f1, 0, 2*pi, 0, 2*pi);
            m(2,1) = integral2(f2, 0, 2*pi, 0, 2*pi);
        end
        
        function pm = angularProductMomentNumerical(this, n)
            % Calculate E(e^inx * e^iny) as used by Jammalamadaka 1988
            %
            % Parameters:
            %   n (scalar)
            %       number of moment
            % Returns:
            %   m (scalar)
            %       n-th angular product moment (complex number)       
            m = this.circularMean();
            f1 = @(x,y) reshape(this.pdf([x(:)'+m(1);y(:)'+m(2)]).*exp(1i*n*x(:)').*exp(1i*n*y(:)'), size(x,1),size(x,2));
            pm = integral2(f1, 0, 2*pi, 0, 2*pi);
        end
        
        function m = circularMean(this)
            % Calculate the circular mean
            %
            % Returns:
            %   m (2 x 1)
            %       circular mean in [0, 2pi)^2
            a = this.trigonometricMoment(1);
            m = mod(atan2(imag(a),real(a)), 2*pi);
        end
        
        function rhoc = circularCorrelationJupp(this)
            % Get the circular correlation coefficient as defined by Jupp
            %
            % Returns:
            %   rhoc (scalar)
            %       Jupp's correlation coefficient
            %
            % see "A General Correlation Coefficient for Directional Data and 
            % Related Regression Problems", P. E. Jupp and K. V. Mardia,
            % Biometrika, Vol. 67, No. 1 (Apr., 1980), pp. 163-173
            C = this.covariance4D();
            % use signed nonsquared version
            rhoc = sign(det(C(1:2,3:4))) * sqrt(trace (C(1:2,1:2) \ C(1:2,3:4) / C(3:4,3:4) * C(3:4,1:2) ));
        end
                
        function rhoc = circularCorrelationJohnson(this)
            % Get the circular correlation coefficient as defined by Johnson
            %
            % Returns:
            %   rhoc (scalar)
            %       Johnson's correlation coefficient
            %
            % see Johnson, R. A. & Wehrly, T. "Measures and Models for Angular 
            % Correlation and Angular-Linear Correlation", 
            % Journal of the Royal Statistical Society. Series B (Methodological), 
            % Wiley for the Royal Statistical Society, 1977, 39, 222-229
            C = this.covariance4D();
            % use signed nonsquared version
            rhoc = sign(det(C(1:2,3:4))) * sqrt(max(eig(C(1:2,1:2)\C(1:2,3:4)/C(3:4,3:4) * C(3:4,1:2) )));
        end        
        
        function mu = mean4D(this)
            % Calculates 4D mean of [cos(x1), sin(x1), cos(x2), sin(x2)]   
            %
            % Returns:
            %   mu (4 x 1)
            %       expectation value of [cos(x1), sin(x1), cos(x2), sin(x2)]
            m = this.trigonometricMoment(1);
            mu = [real(m(1)); imag(m(1)); real(m(2)); imag(m(2))];
        end
        
        function C = covariance4D(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), cos(x2), sin(x2)] 
            %
            % Returns:
            %   C (4 x 4)
            %       covariance matrix of [cos(x1), sin(x1), cos(x2), sin(x2)] 
            C = this.covariance4DNumerical; % fall back to numerical solution
        end
        
        function C = covariance4DNumerical(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), cos(x2), sin(x2)] 
            %
            % Returns:
            %   C (4 x 4)
            %       covariance matrix of [cos(x1), sin(x1), cos(x2), sin(x2)]             
            m = this.mean4D();
            
            f1 = @(x,y) cos(x(:)') - m(1);
            f2 = @(x,y) sin(x(:)') - m(2);
            f3 = @(x,y) cos(y(:)') - m(3);
            f4 = @(x,y) sin(y(:)') - m(4);
            
            f11 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f1(x,y).^2, size(x,1),size(x,2));
            f12 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f1(x,y).*f2(x,y), size(x,1),size(x,2));
            f13 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f1(x,y).*f3(x,y), size(x,1),size(x,2));
            f14 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f1(x,y).*f4(x,y), size(x,1),size(x,2));
            f22 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f2(x,y).^2, size(x,1),size(x,2));
            f23 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f2(x,y).*f3(x,y), size(x,1),size(x,2));
            f24 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f2(x,y).*f4(x,y), size(x,1),size(x,2));
            f33 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f3(x,y).^2, size(x,1),size(x,2));
            f34 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f3(x,y).*f4(x,y), size(x,1),size(x,2));
            f44 = @(x,y) reshape(this.pdf([x(:)';y(:)']).*f4(x,y).^2, size(x,1),size(x,2));
            
            c11 = integral2(f11, 0, 2*pi, 0, 2*pi);
            c12 = integral2(f12, 0, 2*pi, 0, 2*pi);
            c13 = integral2(f13, 0, 2*pi, 0, 2*pi);
            c14 = integral2(f14, 0, 2*pi, 0, 2*pi);
            c22 = integral2(f22, 0, 2*pi, 0, 2*pi);
            c23 = integral2(f23, 0, 2*pi, 0, 2*pi);
            c24 = integral2(f24, 0, 2*pi, 0, 2*pi);
            c33 = integral2(f33, 0, 2*pi, 0, 2*pi);
            c34 = integral2(f34, 0, 2*pi, 0, 2*pi);
            c44 = integral2(f44, 0, 2*pi, 0, 2*pi);
            
            C = [ c11 c12 c13 c14; c12 c22 c23 c24; c13 c23 c33 c34; c14 c24 c34 c44];
        end        
        
        function kld = kldNumerical(this, other)
            % Calculates the Kullback-Leibler divergence between this distribution and another distribution on the Torus
            % Be aware that the kld is not symmetric.
            %
            % Parameters:
            %   other (AbstractCircularDistribution)
            %       distribution to compare to
            % Returns:
            %   kld (scalar)
            %       kld of this distribution to other distribution
            assert(isa(other,'AbstractToroidalDistribution'));
            g = @(x,y) x.*log(x./y);
            f = @(x,y) reshape(g(this.pdf([x(:)';y(:)']),other.pdf([x(:)';y(:)'])), size(x,1),size(x,2));
            kld = integral2(f, 0, 2*pi, 0, 2*pi);
        end
        
        function rhoc = circularCorrelationJammalamadaka(this)
            % Calculates Jammalamadaka's correlation coefficient numerically
            % see Jammalamadaka, S. R. & Sarma, Y. 
            % A Correlation Coefficient for Angular Variables 
            % Statistical Theory and Data Analysis II, 1988, 349-364
            %
            % Returns:
            %   rhoc (scalar)
            %       Jammalamadaka's correlation coefficient            
            rhoc = this.circularCorrelationJammalamadakaNumerical(); % fall back to numerical solution by default
        end
        
        function rhoc = circularCorrelationJammalamadakaNumerical(this)
            % Calculates Jammalamadaka's correlation coefficient numerically
            % see Jammalamadaka, S. R. & Sarma, Y. 
            % A Correlation Coefficient for Angular Variables 
            % Statistical Theory and Data Analysis II, 1988, 349-364
            %
            % Returns:
            %   rhoc (scalar)
            %       Jammalamadaka's correlation coefficient            
            m = this.circularMean();
            fsinAsinB = @(x,y) reshape(this.pdf([x(:)';y(:)']).*sin(x(:)'-m(1)).*sin(y(:)'-m(2)), size(x,1),size(x,2));
            fsinAsquared= @(x,y) reshape(this.pdf([x(:)';y(:)']).*sin(x(:)'-m(1)).^2, size(x,1),size(x,2));
            fsinBsquared= @(x,y) reshape(this.pdf([x(:)';y(:)']).*sin(y(:)'-m(2)).^2, size(x,1),size(x,2));
            EsinAsinB = integral2(fsinAsinB, 0, 2*pi, 0, 2*pi);
            EsinAsquared = integral2(fsinAsquared, 0, 2*pi, 0, 2*pi);
            EsinBsquared = integral2(fsinBsquared, 0, 2*pi, 0, 2*pi);
            rhoc = EsinAsinB/sqrt(EsinAsquared * EsinBsquared);
        end

        function r = integral(this, l1, r1, l2, r2)
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
            r = this.integralNumerical(l1, r1, l2, r2);
        end
        
        function r = integralNumerical(this, l1, r1, l2, r2)
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
            f = @(x,y) reshape(this.pdf([x(:)';y(:)']), size(x,1),size(x,2));
            r = integral2(f, l1, r1, l2, r2);
        end
        
        function l = logLikelihood(this, samples)
            % Calculates the log-likelihood of the given samples
            %
            % Parameters:
            %   sampls (2 x n)
            %       n samples on the torus
            % Returns:
            %   l (scalar)
            %       log-likelihood of obtaining the given samples
            assert(size(samples,1)==2);
            assert(size(samples,2)>=1);
            
            l = sum(log(this.pdf(samples)));
        end   
        
        function s = sample(this, n)
            % Obtain n samples from the distribution
            % use metropolics hastings by default
            % individual distributions can override this with more
            % sophisticated solutions
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (2 x n)
            %       one sample per column
            s = sampleMetropolisHastings(this, n);
        end
                
        function s = sampleMetropolisHastings(this, n)
            % Metropolis Hastings sampling algorithm
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (2 x n)
            %       one sample per column
            %
            % Hastings, W. K. 
            % Monte Carlo Sampling Methods Using Markov Chains and Their Applications 
            % Biometrika, 1970, 57, 97-109
            
            burnin = 10;
            skipping = 5;
            
            totalSamples = burnin+n*skipping;
            s = zeros(2,totalSamples);
            x = this.circularMean; % start with mean
            % A better proposal distribution could be obtained by roughly estimating
            % the uncertainty of the true distribution.
            proposal = @(x) mod(x + mvnrnd([0;0],eye(2,2))', 2*pi); 
            i=1;
            while i<=totalSamples
                xNew = proposal(x); %generate new sample
                a = this.pdf(xNew)/this.pdf(x);
                if a>1
                    %keep sample
                    s(:,i)=xNew;
                    x = xNew;
                    i=i+1;
                else
                    r = rand(1);
                    if a > r
                        %keep sample
                        s(:,i)=xNew;
                        x = xNew;
                        i=i+1;
                    else
                        %reject sample
                    end
                end
            end
            s = s(:,burnin+1:skipping:end);
        end
    end
    
end

