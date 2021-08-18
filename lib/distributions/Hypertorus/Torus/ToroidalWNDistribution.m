classdef ToroidalWNDistribution < AbstractToroidalDistribution & HypertoroidalWNDistribution
    % Bivariate wrapped normal distribution on the torus.
    %
    % see Gerhard Kurz, Igor Gilitschenski, Maxim Dolgov, Uwe D. Hanebeck,
    % Bivariate Angular Estimation Under Consideration of Dependencies Using Directional Statistics
    % Proceedings of the 53rd IEEE Conference on Decision and Control (CDC 2014), Los Angeles, California, USA, December 2014.
    
    methods
        function this = ToroidalWNDistribution (mu_, C_)
            % Constructor
            %
            % Parameters:
            %   mu_ (2 x 1)
            %       mean vector of Gauss before wrapping
            %   C_ (2 x 2)
            %       covariance matrix of Gauss before wrapping
            assert(size(mu_,1)==2,'mu must be 2 x 1'); % Rest of the conditions are tested in superclass constructor
            this@HypertoroidalWNDistribution(mu_, C_);
        end
        
        function p = pdf(this, xa, n)
            % Evaluate pdf at each column of xa.
            %
            % Parameters:
            %   xa (2 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            assert(size(xa, 1) == 2);
            
            %todo: choose n automatically
            if nargin < 3
                n=3;
            end
            
            % call c implementation
            p = toroidalwnpdf(xa,n,this.mu,this.C);
        end
        
        function mu = mean4D(this)
            % Calculates 4D mean of [cos(x1), sin(x1), cos(x2), sin(x2)]   
            %
            % Returns:
            %   mu (4 x 1)
            %       expectation value of [cos(x1), sin(x1), cos(x2), sin(x2)]
            s = this.mu;
            mu = [cos(s(1,:))*exp(-this.C(1,1)/2);
                sin(s(1,:))*exp(-this.C(1,1)/2);
                cos(s(2,:))*exp(-this.C(2,2)/2);
                sin(s(2,:))*exp(-this.C(2,2)/2);
                ];
        end
        
        function C_ = covariance4D(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), cos(x2), sin(x2)]
            %
            % Returns:
            %   C (4 x 4)
            %       covariance matrix of [cos(x1), sin(x1), cos(x2), sin(x2)] 
            %
            % Kurz, Gerhard
            % Directional Estimation for Robotic Beating Heart Surgery
            % Karlsruhe Institute of Technology, 
            % Intelligent Sensor-Actuator-Systems Laboratory, 2015 
            
            C_ (1,1) = 1/2*(1-exp(-this.C(1,1)))*(1-exp(-this.C(1,1))*cos(2*this.mu(1)));
            C_ (1,2) = -1/2*(1-exp(-this.C(1,1)))*exp(-this.C(1,1))*sin(2*this.mu(1));
            C_ (2,1) = C_(1,2);
            C_ (2,2) = 1/2*(1-exp(-this.C(1,1)))*(1+exp(-this.C(1,1))*cos(2*this.mu(1)));
            
            C_ (3,3) = 1/2*(1-exp(-this.C(2,2)))*(1-exp(-this.C(2,2))*cos(2*this.mu(2)));
            C_ (3,4) = -1/2*(1-exp(-this.C(2,2)))*exp(-this.C(2,2))*sin(2*this.mu(2));
            C_ (4,3) = C_(3,4);
            C_ (4,4) = 1/2*(1-exp(-this.C(2,2)))*(1+exp(-this.C(2,2))*cos(2*this.mu(2)));
            
            C_ (1,3) = 1/2*exp(-this.C(1,1)/2-this.C(2,2)/2) * ( exp(-this.C(1,2))*cos(this.mu(1)+this.mu(2)) + exp(this.C(1,2)) * cos(this.mu(1)-this.mu(2)) - 2*cos(this.mu(1))*cos(this.mu(2)));
            C_ (1,4) = 1/2*exp(-this.C(1,1)/2-this.C(2,2)/2) * ( exp(-this.C(1,2))*sin(this.mu(1)+this.mu(2)) - exp(this.C(1,2)) * sin(this.mu(1)-this.mu(2)) - 2*cos(this.mu(1))*sin(this.mu(2)));
            C_ (2,3) = 1/2*exp(-this.C(1,1)/2-this.C(2,2)/2) * ( exp(-this.C(1,2))*sin(this.mu(1)+this.mu(2)) + exp(this.C(1,2)) * sin(this.mu(1)-this.mu(2)) - 2*sin(this.mu(1))*cos(this.mu(2)));
            C_ (2,4) = -1/2*exp(-this.C(1,1)/2-this.C(2,2)/2) * ( exp(-this.C(1,2))*cos(this.mu(1)+this.mu(2)) - exp(this.C(1,2)) * cos(this.mu(1)-this.mu(2)) + 2*sin(this.mu(1))*sin(this.mu(2)));
            
            C_(3,1) = C_(1,3);
            C_(3,2) = C_(2,3);
            C_(4,1) = C_(1,4);
            C_(4,2) = C_(2,4);
        end             
                
        function twn = multiplyMomentBased(this, twn2)
            % Multply two WN distributions (approximate using moment
            % matching).
            %
            % Parameters:
            %   twn2 (ToroidalWNDistribution)
            %       distribution to multiply with
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       product of this and twn2
            
            % get moments from 1D WN  using numerical integration
            f0 = @(x,y) reshape( this.pdf([x(:)';y(:)']) .* twn2.pdf([x(:)';y(:)']) , size(x,1),size(x,2));
            c=integral2(f0, 0, 2*pi, 0, 2*pi);
            f1 = @(x,y) reshape( exp(1i*x(:)').* this.pdf([x(:)';y(:)']) .* twn2.pdf([x(:)';y(:)']) , size(x,1),size(x,2));
            m1 = integral2(f1, 0, 2*pi, 0, 2*pi)/c;
            f2 = @(x,y) reshape( exp(1i*y(:)').* this.pdf([x(:)';y(:)']) .* twn2.pdf([x(:)';y(:)']) , size(x,1),size(x,2));
            m2 = integral2(f2, 0, 2*pi, 0, 2*pi)/c;
            mu_ = mod([atan2(imag(m1), real(m1)); atan2(imag(m2), real(m2))],2*pi);
            si1squared = -2 * log(abs(m1));
            si2squared = -2 * log(abs(m2));
            
            % get correlation
            f1 = @(x,y) reshape( sin(x(:)'-mu_(1)).*sin(y(:)'-mu_(2)).* this.pdf([x(:)';y(:)']) .* twn2.pdf([x(:)';y(:)']) , size(x,1),size(x,2));
            a = integral2(f1, 0, 2*pi, 0, 2*pi);
            f2 = @(x,y) reshape( sin(x(:)'-mu_(1)).^2 .* this.pdf([x(:)';y(:)']) .* twn2.pdf([x(:)';y(:)']) , size(x,1),size(x,2));
            b = integral2(f2, 0, 2*pi, 0, 2*pi);
            f3 = @(x,y) reshape( sin(y(:)'-mu_(2)).^2 .* this.pdf([x(:)';y(:)']) .* twn2.pdf([x(:)';y(:)']) , size(x,1),size(x,2));
            c = integral2(f3, 0, 2*pi, 0, 2*pi);
            circCorrelation = a/sqrt(b*c); %a/sqrt(b*c) is the circular correlation coefficient
            z = sqrt(sinh(si1squared)*sinh(si2squared))* circCorrelation;
            %function for inverse of sinh
            arsinh = @(z) log(z + sqrt(z^2+1));
            si12 = arsinh(z);
            
            rho = si12/sqrt(si1squared*si2squared);
            if abs(rho)>=1
                error('conversion impossible')
            end
            
            % return result
            C_ = [si1squared si12;
                si12 si2squared];
            twn = ToroidalWNDistribution(mu_, C_);
        end
        
        function twn = convolve(this, twn2)
            % Calculate convolution of two ToroidalWNDistribution objects
            % (exact).
            %
            % Parameters:
            %   twn2 (ToroidalWNDistribution)
            %       distribution to convolve with
            % Returns:
            %   twn (ToroidalWNDistribution)
            %       convolution of this and twn2
            
            assert(isa(twn2, 'ToroidalWNDistribution'));
            
            mu_ = mod(this.mu+twn2.mu,2*pi);
            C_ = this.C + twn2.C;
            twn = ToroidalWNDistribution(mu_, C_);
        end

        function rho = realCorrelation(this)
            % Get the correlation coefficient from the covariance matrix.
            %
            % Returns:
            %   rho (scalar)
            %       Pearson correlation coefficient of a Gaussian with identical
            %       parameters
            rho = this.C(1,2)/sqrt(this.C(1,1)*this.C(2,2));
        end
        
        function rhoc = circularCorrelationJammalamadaka(this)
            % Get the circular correlation coefficient as defined by
            % Jammalamadaka.
            %
            % Jammalamadaka, S. R. & Sarma, Y. A
            % Correlation Coefficient for Angular Variables
            % Statistical Theory and Data Analysis II, 1988, 349-364
            %
            % Returns:
            %   rhoc (scalar)
            %       Jammalamadaka's correlation coefficient                      
            rhoc = sinh(this.C(1,2))/sqrt(sinh(this.C(1,1))*sinh(this.C(2,2)));
        end
                                
        function tvm = toToroidalVMSine(this)
            % Convert to Torus VM (sine version).
            %
            % Returns:
            %   tvm (ToroidalVMSineDistribution)
            %       ToroidalVMSineDistribution that approximates this
            %       distribution
            %
            % obtained by inverting formulas from 
            % Singh, H.; Hnizdo, V. & Demchuk, E. 
            % Probabilistic Model for Two Dependent Circular Variables 
            % Biometrika, 2002, 89, 719-723
            si1squared = this.C(1,1);
            si2squared = this.C(2,2);
            rho = this.C(1,2)/sqrt(si1squared*si2squared);
            kappa(1,1) = 1/si1squared/(1-rho^2);
            kappa(2,1) = 1/si2squared/(1-rho^2);
            lambda = rho/sqrt(si1squared*si2squared*(rho^2-1)^2);
            tvm = ToroidalVMSineDistribution(this.mu, kappa, lambda);
        end
        
        function wn = marginalizeTo1D(this, dimension)
            % Get marginal distribution in first or second dimension, i.e., 
            % f(x_1) or f(x_2), respectively.
            %
            % Parameters:
            %   dimension (scalar)
            %       the marginal in which dimension to calculate (1 or 2),
            %       the other dimension is marginalized out
            % Returns:
            %   wn (WNDistribution)
            %       marginal distribution (marginals are WN-distributed)
            assert(dimension==1 || dimension==2);
            wn = WNDistribution(this.mu(dimension), sqrt(this.C(dimension, dimension)));
        end
    end
    
    methods (Static)
        function twn = mleJensen(samples)
            % Obtain WN distribution from samples with approximate MLE using 
            % Jensen's inqueality.
            %
            % Parameters:
            %   samples (2 x n matrix)
            %       samples on the torus
            % Returns: 
            %   twn (ToroidalWNDistribution)
            %       distribution obtained with approximate MLE
            
            % based on
            % Expectation Maximization - Math and Pictures,
            % Johannes Traa, 2014
            assert(size(samples,1)==2);
            assert(size(samples,2)>=3);
            terms = 1;
                                   
            %start value
            twd = ToroidalWDDistribution(samples);
            mualpha = twd.circularMean();
            xbar = sum( twd.w.* cos(twd.d(1,:) - mualpha(1)));
            ybar = sum( twd.w.* cos(twd.d(2,:) - mualpha(2)));
            c11 = -2 * log(xbar);
            c22 = -2 * log(ybar);
            c12 = 0; 
            sigmaalpha = [c11, c12; c12, c22];
            
            function a = alphaNumerator(w1, w2, s)
                %a = mvnpdf(s+[2*w1*pi; 2*w2*pi], mualpha, sigmaalpha);
                a = mvnpdffast(s'+[2*w1*pi, 2*w2*pi], mualpha', sigmaalpha);
                %a = mvnpdffast_mex(s'+[2*w1*pi, 2*w2*pi], mualpha', sigmaalpha);
            end
            
            n = size(samples,2);
            mu_ = 0;
            C_ = 0;
            alphaDenominator = zeros(1,n);
            for i=1:n
                alphaDenominator(i) = toroidalwnpdf(samples(:,i),terms,mualpha,sigmaalpha);
                for w1=-terms:terms
                    for w2=-terms:terms
                        mu_ = mu_ + alphaNumerator(w1,w2,samples(:,i))/alphaDenominator(i)*(samples(:,i)+[2*w1*pi; 2*w2*pi]);
                    end
                end
            end
            mu_=mu_/n;
            for i=1:n
                for w1=-terms:terms
                    for w2=-terms:terms
                        C_ = C_ + alphaNumerator(w1, w2, samples(:,i))/alphaDenominator(i)*(samples(:,i)+[2*w1*pi; 2*w2*pi] - mu_)*(samples(:,i)+[2*w1*pi; 2*w2*pi] - mu_)';
                    end
                end
            end
            C_ = C_/n;
            C_ = (C_+ C_')/2; %ensure symmetry in presence of numerical errors
            twn = ToroidalWNDistribution(mu_,C_);
        end
        
        function twn = mleNumerical(samples)
            % Obtain TWN distribution from samples using MLE (based on
            % numerical optimization).
            %
            % Parameters:
            %   samples (2 x n matrix)
            %       samples on the torus
            % Returns: 
            %   twn (ToroidalWNDistribution)
            %       distribution obtained with numerical MLE            
            function y = f(x)
                Ctemp = [x(3) x(4); x(4) x(5)];
                if x(3)>0 && det(Ctemp)>0
                    w = ToroidalWNDistribution([x(1);x(2)],Ctemp);
                    y = -w.logLikelihood(samples);
                else
                    y = Inf;
                end
            end
            
            %start value
            twd = ToroidalWDDistribution(samples);
            mu_ = twd.circularMean();
            xbar = sum( twd.w.* cos(twd.d(1,:) - mu_(1)));
            ybar = sum( twd.w.* cos(twd.d(2,:) - mu_(2)));
            c11 = -2 * log(xbar);
            c22 = -2 * log(ybar);
            c12 = 0; 
            startValue = [mu_(1), mu_(2), c11, c12, c22];
            
            %perform optimization
            x = fminsearch(@f,startValue,optimset('display','none'));
            twn = ToroidalWNDistribution([x(1);x(2)],[x(3) x(4); x(4) x(5)]);
        end
    end
end
