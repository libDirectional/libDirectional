classdef SE2PWNDistribution
    % Partially wrapped normal distribution for SE(2). The first component 
    % is wrapped, the second and the third are not.
    % 
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % The Partially Wrapped Normal Distribution for SE(2) Estimation
    % Proceedings of the 2014 IEEE International Conference on Multisensor
    % Fusion and Information Integration (MFI 2014), Beijing, China, September 2014.
    %
    % Gerhard Kurz
    % Directional Estimation for Robotic Beating Heart Surgery
    % Karlsruhe Institute of Technology, 
    % Intelligent Sensor-Actuator-Systems Laboratory, 2015 

    properties
        mu
        C
    end
    
    methods
        function this = SE2PWNDistribution(mu_, C_)
            % Constructor
            % check parameters
            assert(all(size(mu_) == [3,1]), 'mu must 3x1'); 
            assert(all(size(C_) == [3,3]), 'C must be 3x3');
            assert(all(all(C_ == C_')), 'C must be symmetric');
            assert(all(eig(C_) > 0), 'C must be positive definite');
            
            % assign parameters
            this.mu = mu_;
            this.mu(1) = mod(this.mu(1),2*pi);
            this.C= C_;
        end
        
        function p = pdf(this, xa)
            % Evaluate pdf at each column of xa
            %
            % Parameters:
            %   xa (3 x n)
            %       n locations where to evaluate the pdf
            % Returns:
            %   p (1 x n)
            %       value of the pdf at each location
            
            assert(size(xa,1)==3);
            
            if size(xa,2)>1
                p = zeros(1,size(xa,2));
                for i=1:size(xa,2)
                    p(i) = this.pdf(xa(:,i));
                end
                return
            end
            
            Cinv = -0.5 * inv(this.C); %#ok<MINV>
            xa = xa-this.mu;
            n = 1;
            sum = exp(xa'*Cinv*xa); %0 th summand
            while true
                dxa = [2*pi*n;0; 0];
                sumNew = sum + exp((xa+dxa)'*Cinv*(xa+dxa)) + exp((xa-dxa)'*Cinv*(xa-dxa)); %-nth and +nth summand
                if sumNew == sum %break when result is accurate according to machine precision
                    break;
                end
                sum = sumNew;
                n = n + 1;
            end
            p = sum/( (2*pi)^(3/2) * sqrt(det(this.C)) );
        end
        
        function plot3d(this)
            % Creates a three-dimensional plot
            stepCirc = 0.5;
            stepLin = 0.5;
            linSize = 3; 
            x2min = this.mu(2)-linSize*sqrt(this.C(2,2));
            x2max = this.mu(2)+linSize*sqrt(this.C(2,2));
            x3min = this.mu(3)-linSize*sqrt(this.C(3,3));
            x3max = this.mu(3)+linSize*sqrt(this.C(3,3));
            pdfMax = this.pdf(this.mu);
            [X,Y,Z] = sphere(4);
            clf 
            hold on
            for x1=0:stepCirc:2*pi
                for x2=x2min:stepLin:x2max
                    for x3=x3min:stepLin:x3max
                        ratio = this.pdf([x1; x2; x3])/pdfMax;
                        plotsize = 0.5*ratio;
                        color = jet;
                        color = color(1+floor(ratio*63),:);
                        if plotsize>0.01
                            surf(plotsize*X+x1,plotsize*Y+x2,plotsize*Z+x3, 'facecolor', color);
                        end
                    end
                end
            end
            hold off
            xlabel('x_1')
            setupAxisCircular('x')
            ylabel('x_2')
            zlabel('x_3')
            view(40,20)
            grid
        end
                                     
        function pwn = convolve(this, pwn2)
            % Calculate convolution of two PWN distributions (exact)
            %
            % Parameters:
            %   pwn2 (SE2PWNDistribution)
            %       distribution to convolve with
            % Returns:
            %   pwn (SE2PWNDistribution)
            %       convolution of this and pwn2
            
            assert(isa(pwn2, 'SE2PWNDistribution'));
            
            mu_ = this.mu+pwn2.mu;
            mu_(1) = mod(mu_(1),2*pi);
            C_ = this.C + pwn2.C;
            pwn = SE2PWNDistribution(mu_, C_);
        end
        
        function mu = mean4D(this)
            % Calculates 4D mean of [cos(x1), sin(x1), x2, x3]
            %
            % Returns:
            %   mu (4 x 1)
            %       expectation value of [cos(x1), sin(x1), x2, x3]
            s = this.mu;
            mu = [cos(s(1,:))*exp(-this.C(1,1)/2); sin(s(1,:))*exp(-this.C(1,1)/2); s(2:3,:)];
        end
        
        function C = covariance4D(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), x2, x3]
            %
            % Returns:
            %   C (4 x 4)
            %       covariance matrix of [cos(x1), sin(x1), x2, x3]
            %
            % Generalization of Johnson 1977, page 224
            %
            % Kurz, Gerhard
            % Directional Estimation for Robotic Beating Heart Surgery
            % Karlsruhe Institute of Technology, 
            % Intelligent Sensor-Actuator-Systems Laboratory, 2015 

            a = 1 - exp(-this.C(1,1));
            si11 = 1/2*a*(1-exp(-this.C(1,1))*cos(2*this.mu(1)));
            si22 = 1/2*a*(1+exp(-this.C(1,1))*cos(2*this.mu(1)));
            si33 = this.C(2,2);
            si44 = this.C(3,3);
            
            si12 = -a/2*exp(-this.C(1,1))*sin(2*this.mu(1));
            
            si13 = -exp(-1/2*this.C(1,1))*this.C(1,2)*sin(this.mu(1));
            si23 = exp(-1/2*this.C(1,1))*this.C(1,2)*cos(this.mu(1));
            si14 = -exp(-1/2*this.C(1,1))*this.C(1,3)*sin(this.mu(1));
            si24 = exp(-1/2*this.C(1,1))*this.C(1,3)*cos(this.mu(1));
            si34 = this.C(2,3);
            
            C = [si11 si12 si13 si14; 
                 si12 si22 si23 si24;
                 si13 si23 si33 si34;
                 si14 si24 si34 si44];
        end
        
        function C = covariance4DNumerical(this)
            % Calculates 4D covariance of [cos(x1), sin(x1), x2, x3]
            % numerically
            %
            % Returns:
            %   C (4 x 4 matrix)
            %       covariance matrix of [cos(x1), sin(x1), x2, x3]
            s = this.sample(1000);
            S = [cos(s(1,:)); sin(s(1,:)); s(2:3,:)];
            C = cov(S');
            
            % An alternative solution based on numerical integration would
            % be possible, similar to the toroidal case.
        end
        
        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (3 x n matrix)
            %       n samples on [0,2pi) x R^2
            %
            % sample multivariate normal distribution, then wrap
            assert(isscalar(n));
            s = mvnrnd(this.mu, this.C, n)';
            s(1,:) = mod(s(1,:),2*pi);
        end

        function gauss = toGaussian(this)
            % Convert to 3D Gaussian
            %
            % Returns:
            %   gauss (GaussianDistribution)
            %       GaussianDistribution with same parameters
            %
            % this is a simple conversion that just keeps the parameters
            % for large uncertanties, better conversions are possible
            gauss = GaussianDistribution(this.mu, this.C);
        end
               
        function gauss = marginalizeCircular(this)
            gauss = GaussianDistribution(this.mu(2:3), this.C(2:3,2:3));
        end
        
        function wn = marginalizeLinear(this)
            wn = WNDistribution(this.mu(1), sqrt(this.C(1,1)));
        end
        
        function r = integral(this)
            % Integrates the density to check normalization
            %
            % Returns:
            %   r (scalar)
            %       integral over the entire density
            %
            
            % integrate circular part from 0 to 2pi
            l1 = 0;
            r1 = 2*pi;
            % The linear part should be integrated from -Inf to Inf, but
            % Matlab's numerical integration does not handle that well.
            % For this reason we integrate from mu-10*sigma to mu+10*sigma,
            % which contains almost the entire probability mass
            intervalSize = 10;
            si2 = sqrt(this.C(2,2));
            l2 = this.mu(2)-intervalSize*si2;
            r2 = this.mu(2)+intervalSize*si2;
            si3 = sqrt(this.C(3,3));
            l3 = this.mu(3)-intervalSize*si3;
            r3 = this.mu(3)+intervalSize*si3;
            
            f = @(x,y,z) reshape(this.pdf([x(:)';y(:)';z(:)']), size(x,1), size(x,2));
            r = integral3(f, l1, r1, l2, r2, l3, r3);
        end      
    end
    
    methods (Static)
        function pwn = fromSamples(samples)
            % Estimate pwn from samples using moments
            %
            % Parameters:
            %   samples (3 x n matrix)
            %       samples on [0,2pi) x R^2
            % Returns:
            %   pwn (SE2PWNDistribution)
            %       fitted PWN distribution
            assert(size(samples,1)==3);
            
            S = [cos(samples(1,:)); sin(samples(1,:)); samples(2:3,:)]; %enhance to 4D
            mu4 = mean(S,2);
            mu_(1) = mod(atan2(mu4(2),mu4(1)),2*pi); %argument of first circular moment
            m1abs = sqrt(mu4(1)^2 + mu4(2)^2); %absolute value of first circular moment
            mu_(2:3) = mu4(3:4);
            
            C4 = cov(S'); %4D covariance
            
            C_ = zeros(3,3);
            C_(1,1) = -2 * log(m1abs); %sigma^2 of WN
            
            %circular linear correlation (minimizes squred error)
            C_(1,2) = -C4(1,3)*exp(1/2*C_(1,1))*sin(mu_(1)) + C4(2,3)*exp(1/2*C_(1,1))*cos(mu_(1));
            C_(1,3) = -C4(1,4)*exp(1/2*C_(1,1))*sin(mu_(1)) + C4(2,4)*exp(1/2*C_(1,1))*cos(mu_(1)); 
            
            %mirror
            C_(2,1) = C_(1,2);
            C_(3,1) = C_(1,3);
            
            %copy Gaussian part
            C_(2:3,2:3) = C4(3:4,3:4);
            
            %construct new PWN
            pwn = SE2PWNDistribution(mu_',C_);
        end
    end
    
end
