classdef VMFDistribution < AbstractHypersphericalDistribution
    % Represents a von Mises-Fisher distribution.
    % Notation:
    % In this class, d represents the dimension of the distribution.
    %
    % see R. Fisher, "Dispersion on a sphere," Proceedings of the Royal Society of
    % London. Series A, Mathematical and Physical Sciences, vol. 217, no. 1130,
    % pp. 295-305, 1953.
    
    properties
        kappa (1,1) double  % concentration (scalar)
        mu (:,1) double      % mean as a unit vector
        C {mustBeNonnegative}    % normalization constant
    end
    
    methods
        function VMF = VMFDistribution(mu_, kappa_)
            % Constructor
            %
            % Parameters:
            %   mu_ (d x 1)
            %       location parameter (unit vector)
            %   kappa_ (scalar)
            %       concentration parameter (>=0)
            arguments
                mu_ (:,1) double
                kappa_ (1,1) double {mustBeNonnegative}
            end
            epsilon = 1E-6;
            assert(size(mu_,1) >= 2, 'mu must be at least two-dimensinal for the circular case');
            assert(abs(norm(mu_) - 1)<epsilon, 'mu must be a normalized');
            
            VMF.mu = mu_;
            VMF.kappa = kappa_;
            
            VMF.dim = size(mu_,1);
            if VMF.dim==3
                VMF.C = kappa_/(4*pi*sinh(kappa_));
            else
                VMF.C = kappa_^(VMF.dim/2-1)/ ( (2*pi)^(VMF.dim/2) * besseli(VMF.dim/2-1, kappa_));
            end
        end
        
        function p = pdf(this, xa)
            % Evaluates pdf at each column of xa
            % Parameters:
            %   xa (d x n matrix)
            %       each column represents one of the n points in R^d that the
            %       pdf is evaluated at; can be just a (d x 1) vector as well
            % Returns:
            %   p (1 x n row vector)
            %       values of the pdf at each column of xa
            arguments
                this (1,1) VMFDistribution
                xa (:,:) double
            end
            assert (size(xa,1) == size(this.mu,1));
            
            p = this.C * exp( this.kappa * this.mu' * xa);
        end
        
        function mu = meanDirection(this)
            arguments
                this (1,1) VMFDistribution
            end
           mu = this.mu;
        end
        
        function samples = sampleDeterministic(this)
            % Deterministic sampling based on moment matching
            % Returns:
            %   samples (d x n matrix)
            %       one sample per column
            %
            % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
            % Unscented von Mises-Fisher Filtering
            % IEEE Signal Processing Letters, 2016.
            arguments
                this (1,1) VMFDistribution
            end
            % get samples for vm with mu=[1 0 ... 0] first            
            samples = zeros(this.dim, this.dim*2-1);
            samples(1,1) = 1;
            m1 = besseli(this.dim/2, this.kappa,1)/besseli(this.dim/2-1, this.kappa,1);
            for i=1:this.dim-1
                alpha = acos( (( this.dim*2-1)*m1 - 1)/( this.dim*2-2));
                samples(1,2*i) = cos(alpha);
                samples(1,2*i+1) = cos(alpha);
                samples(i+1,2*i) = sin(alpha);
                samples(i+1,2*i+1) = -sin(alpha);
            end
            
            % rotate samples around mu
            Q = this.getRotationMatrix();
            samples = Q*samples;
        end
        
        function Q = getRotationMatrix(this)
            % Computes a rotation matrix that rotates [1,0,...,0] to mu.
            %
            % Returns:
            %   Q (d x d matrix)
            %       rotation matrix
            %
            % Gerhard Kurz, Uwe D. Hanebeck,
            % Stochastic Sampling of the Hyperspherical von Mises-Fisher Distribution Without Rejection Methods
            % Proceedings of the IEEE ISIF Workshop on Sensor Data Fusion: Trends, Solutions, Applications (SDF 2015), Bonn, Germany, October 2015
            arguments
                this (1,1) VMFDistribution
            end
            M = zeros(this.dim,this.dim);
            M(:,1) = this.mu; 
            [Q,R] = qr(M);
            if R(1,1)<0
                Q=-Q;
            end
            % Q may have det -1, but we ignore that for now because the
            % distribution does not change when mirrored
        end
        
        function m = mode(this)
            % Calculate the mode of a VMF distribution
            % Returns:
            %   m (column vector)
            %       mode of the distribution
            arguments
                this (1,1) VMFDistribution
            end
            m = this.mu;
        end        
        
        function dist = setMode(this, newMode)
            arguments
                this (1,1) VMFDistribution
                newMode (:,1) double
            end
            assert(isequal(size(newMode),size(this.mu)));
            dist = this;
            dist.mu = newMode;
        end

        function vmf = multiply(this, other)
            % Multiplies this density with another VMF density (exact).
            %
            % Parameters:
            %   other (VMFDistribution)
            %       the other VMF distribution
            % Returns:
            %   vmf (VMFDistribution)
            %       the vmf distribution of the renormalized product.
            arguments
                this (1,1) VMFDistribution
                other (1,1) VMFDistribution
            end
            assert(isequal(size(this.mu),size(other.mu)));
            
            mu_ = this.kappa*this.mu + other.kappa*other.mu;
            kappa_ = norm(mu_);
            mu_ = mu_ /kappa_;
            vmf = VMFDistribution(mu_,kappa_);
        end
        
        function vmf = convolve(this, other)
            % Computes a convolution-like operation of this density with 
            % another VMF density (approximate).
            %
            % This is not a true convolution, but just a similar
            % operation. This operation is not symmetric as the mu-value of the second
            % density is ignored. To ensure that the user is aware of this,
            % the last entry of mu has to be set to 1 and all others to 0.
            %
            % Traa, J. & Smaragdis, P. 
            % Multiple Speaker Tracking With the Factorial von Mises--Fisher Filter 
            % IEEE International Workshop on Machine Learning for Signal Processing (MLSP), 2014
            % (see Section 4.1)
            %
            % Parameters:
            %   other (VMFDistribution)
            %       the other VMF distribution
            % Returns:
            %   vmf (VMFDistribution)
            %       the vmf distribution of the renormalized product.     
            arguments
                this (1,1) VMFDistribution
                other (1,1) VMFDistribution
            end
            assert(other.mu(end)==1,'Other is not zonal');
            assert(all(size(this.mu) == size(other.mu)));
            d = this.dim;
            
            mu_ = this.mu;
            kappa_ = VMFDistribution.AdInverse(d, VMFDistribution.Ad(d, this.kappa)*VMFDistribution.Ad(d, other.kappa));
            vmf = VMFDistribution(mu_,kappa_);
        end
        
        function vm = toVM(this)
            % Converts the VMF distribution to a VM for d=2 (this
            % distribution is identical to the VMF, but parameterized
            % using an angle rather than a unit vector).
            % 
            % Returns:
            %   vm (VMDistribution)
            %       the corresponding VM distribution
            arguments
                this (1,1) VMFDistribution
            end
            assert(this.dim==2, 'Conversion only possible in two dimensions.')
            mu_ = mod(atan2(this.mu(2), this.mu(1)),2*pi);
            vm = VMDistribution(mu_, this.kappa);
        end
        
        function g = toGaussian(this)
            % Converts the VMF distribution to a Gaussian by moment
            % matching.
            %
            % Returns:
            %   g (GaussianDistribution)
            %       resulting Gaussian
            arguments
                this (1,1) VMFDistribution
            end
            n = 10000;
            s = this.sample(n);
            %C_ = cov(s');
            C_ = (s-repmat(this.mu,1,n))*(s-repmat(this.mu,1,n))'/(n-1);
            g = GaussianDistribution(this.mu, C_);
        end
        
        function s = sample(this, n)
            % Generate samples from VMF distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   X (d x n matrix)
            %       generated samples (one sample per column)
            arguments
                this (1,1) VMFDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            if this.dim == 2
                s = this.toVM().sample(n);
                s = [cos(s);sin(s)];
            elseif this.dim == 3
                s = this.sampleJakob(n);
            elseif this.dim == 5 || this.dim == 7
                s = this.sampleCdfVectorized(n);
            else
                s = this.sampleWood(n);
            end
        end
        
        function s = sampleJakob(this, n)
            % Generate samples from VMF distribution using Jakob's algorithm
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   X (d x n matrix)
            %       generated samples (one sample per column)
            %
            % Jakob, W. 
            % Numerically Stable Sampling of the von Mises Fisher Distribution on S2 (and other Tricks) 
            % 2015
            %
            % see also
            % Jung, S. 
            % Generating von Mises--Fisher Distribution on the Unit Sphere (S2) 
            % 2009
            arguments
                this (1,1) VMFDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            assert(this.dim==3, 'This method only works in three dimensions.');
            
            Q = this.getRotationMatrix();
            
            U = rand(1,n);
            W = 1 + 1/this.kappa * log(U + (1-U)*exp(-2*this.kappa)); % numerically stable version by Jakob
            % generate V uniform on the unit circle
            alpha = 2*pi*rand(1,n);
            V = [cos(alpha); sin(alpha)];
            % we differ from Jakob's algorithm here and sample a
            % distribution with mu = [1,0,...,0] rather than mu =
            % [0,...,0,1]
            X = [W; sqrt(1-W.^2).*V(1,:); sqrt(1-W.^2).*V(2,:)]; 
            s = Q*X;
        end
        
        function s = sampleWood(this, n)
            % Generate samples from VMF distribution using Wood's algorithm
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   X (d x n matrix)
            %       generated samples (one sample per column)
            %
            % Wood, A. T. A. 
            % Simulation of the von Mises--Fisher distribution 
            % Communications in Statistics---Simulation and Computation, 
            % Taylor & Francis, 1994, 23, 157-164
            %
            % Works for an arbitrary number of dimensions.
            arguments
                this (1,1) VMFDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            lambda = this.kappa;
            m = this.dim;
            s = zeros(m,n);
            Q = this.getRotationMatrix();

            % Step 0
            b = (-2 *lambda + sqrt(4*lambda^2 + (m-1)^2))/(m-1);
            x0 = (1-b)/(1+b);
            c = lambda*x0 + (m-1) * log ( 1-x0^2);

            for i=1:n
                while true
                    % Step 1
                    Z = betarnd((m-1)/2, (m-1)/2); 
                    U = rand(1);
                    W = (1-(1+b)*Z)/(1-(1-b)*Z);

                    % Step 2
                    if lambda*W + (m-1)*log(1-x0*W)-c >= log(U)
                        break
                    end
                end

                % Step 3
                V = mvnrnd(zeros(m-1,1), eye(m-1,m-1));
                V = V'/norm(V);
                % we differ from Wood's algorithm here and sample a
                % distribution with mu = [1,0,...,0] rather than mu =
                % [0,...,0,1]
                X = [W; sqrt(1-W^2)*V]; 

                % rotate samples according to mu
                s(:,i) = Q*X;
            end
        end        
        
        function s = sampleCdf(this, n)
            % Generate samples from VMF distribution using CDF algorithm
            % (non-vectorized version, slower)            
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   X (d x n matrix)
            %       generated samples (one sample per column)
            %
            % Gerhard Kurz, Uwe D. Hanebeck,
            % Stochastic Sampling of the Hyperspherical von Mises-Fisher Distribution Without Rejection Methods
            % Proceedings of the IEEE ISIF Workshop on Sensor Data Fusion: Trends, Solutions, Applications (SDF 2015), Bonn, Germany, October 2015
            arguments
                this (1,1) VMFDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            s = zeros(this.dim,n);
            Q = this.getRotationMatrix();
            k = this.kappa;
            pdf = @(x) exp(this.kappa*cos(x)).*sin(x).^(this.dim-2);

            switch this.dim
                % use 
                % d = 3; syms k x y; int(exp(k*cos(x)).*sin(x).^(d-2), x, 0, y) 
                % to obtain solution for a certain (odd) dimension
                case 3
                    cdf = @(y) -(exp(k*cos(y)) - exp(k))/k;
                case 5
                    cdf = @(y) (exp(k)*(2*k - 2))/k^3 - (exp(k*cos(y))*(k^2*sin(y)^2 + 2*k*cos(y) - 2))/k^3;
                case 7
                    cdf = @(y) (exp(k)*(8*k^2 - 24*k + 24))/k^5 - (exp(k*cos(y))*(8*k^2*cos(y)^2 - 4*k^2*sin(y)^2 + k^4*sin(y)^4 - 24*k*cos(y) + 4*k^3*cos(y)*sin(y)^2 + 24))/k^5;
                case 9  
                    cdf = @(y) (exp(k*cos(y))*(72*k^3*cos(y) - 6*k^5*cos(y) + 360*k^2*cos(y)^2 - 120*k^3*cos(y)^3 - 36*k^4*cos(y)^2 + 30*k^4*cos(y)^4 + 12*k^5*cos(y)^3 + 3*k^6*cos(y)^2 - 6*k^5*cos(y)^5 - 3*k^6*cos(y)^4 + k^6*cos(y)^6 - 720*k*cos(y) - 72*k^2 + 6*k^4 - k^6 + 720))/k^7 + (exp(k)*(48*k^3 - 288*k^2 + 720*k - 720))/k^7;
                case 11
                    cdf = @(y) (384*exp(k))/k^5 - (40320*exp(k*cos(y)))/k^9 - (3840*exp(k))/k^6 + (17280*exp(k))/k^7 - (40320*exp(k))/k^8 + (40320*exp(k))/k^9 + (40320*exp(k*cos(y))*cos(y))/k^8 - (2880*exp(k*cos(y))*(7*cos(y)^2 - 1))/k^7 - (exp(k*cos(y))*(cos(y)^2 - 1)^4)/k - (48*exp(k*cos(y))*(35*cos(y)^4 - 30*cos(y)^2 + 3))/k^5 + (960*exp(k*cos(y))*cos(y)*(7*cos(y)^2 - 3))/k^6 + (8*exp(k*cos(y))*cos(y)*(cos(y)^2 - 1)^3)/k^2 - (8*exp(k*cos(y))*(cos(y)^2 - 1)^2*(7*cos(y)^2 - 1))/k^3 + (48*exp(k*cos(y))*cos(y)*(7*cos(y)^4 - 10*cos(y)^2 + 3))/k^4;
                otherwise
                    warning('falling back to numerical integration')
                    cdf = @(y) integral(pdf, 0, y);
            end
            
            u = rand(1,n)*cdf(pi); % consider the fact that pdf is not normalized
            
            for i=1:n                                                             
                %bisection
                lower = 0;
                upper = pi;
                for j=1:52
                    middle = (lower+upper)/2;
                    cdfMiddle = cdf(middle);
                    if cdfMiddle>u(i)
                        upper = middle;
                    else
                        lower = middle;
                    end
                end
                phi = middle;
                
                
                W = cos(phi);                
                
                % Step 3 from Wood
                V = mvnrnd(zeros(this.dim-1,1), eye(this.dim-1,this.dim-1));
                V = V'/norm(V);
                % we differ from Wood's algorithm here and sample a
                % distribution with mu = [1,0,...,0] rather than mu =
                % [0,...,0,1]
                X = [W; sqrt(1-W^2)*V]; 
                
                s(:,i) = Q*X;
            end
        end
        
        function s = sampleCdfVectorized(this, n)
            % Generate samples from VMF distribution using CDF algorithm
            % (vectorized version, faster)
            %
            % Parameters:
            %   n (scalar)
            %       number of samples to generate
            % Returns:
            %   X (d x n matrix)
            %       generated samples (one sample per column)
            %
            % Gerhard Kurz, Uwe D. Hanebeck,
            % Stochastic Sampling of the Hyperspherical von Mises-Fisher Distribution Without Rejection Methods
            % Proceedings of the IEEE ISIF Workshop on Sensor Data Fusion: Trends, Solutions, Applications (SDF 2015), Bonn, Germany, October 2015
            arguments
                this (1,1) VMFDistribution
                n (1,1) {mustBePositive,mustBeInteger}
            end
            s = zeros(this.dim,n);
            Q = this.getRotationMatrix();
            k = this.kappa;
            pdf = @(x) exp(this.kappa*cos(x)).*sin(x).^(this.dim-2);

            switch this.dim
                % use 
                % d = 3; syms k x y; int(exp(k*cos(x)).*sin(x).^(d-2), x, 0, y) 
                % to obtain solution for a certain (odd) dimension
                case 3
                    cdf = @(y) -(exp(k*cos(y)) - exp(k))/k;
                case 5
                    cdf = @(y) (exp(k)*(2*k - 2))/k^3 - (exp(k*cos(y)).*(k^2*sin(y).^2 + 2*k*cos(y) - 2))/k^3;
                case 7
                    cdf = @(y) (exp(k)*(8*k^2 - 24*k + 24))/k^5 - (exp(k*cos(y)).*(8*k^2*cos(y).^2 - 4*k^2*sin(y).^2 + k^4*sin(y).^4 - 24*k*cos(y) + 4*k^3*cos(y).*sin(y).^2 + 24))/k^5;
                case 9  
                    cdf = @(y) (exp(k*cos(y)).*(72*k^3*cos(y) - 6*k^5*cos(y) + 360*k^2*cos(y).^2 - 120*k^3*cos(y).^3 - 36*k^4*cos(y).^2 + 30*k^4*cos(y).^4 + 12*k^5*cos(y).^3 + 3*k^6*cos(y).^2 - 6*k^5*cos(y).^5 - 3*k^6*cos(y).^4 + k^6*cos(y).^6 - 720*k*cos(y) - 72*k^2 + 6*k^4 - k^6 + 720))/k^7 + (exp(k)*(48*k^3 - 288*k^2 + 720*k - 720))/k^7;
                case 11
                    cdf = @(y) (384*exp(k))/k^5 - (40320*exp(k*cos(y)))/k^9 - (3840*exp(k))/k^6 + (17280*exp(k))/k^7 - (40320*exp(k))/k^8 + (40320*exp(k))/k^9 + (40320*exp(k*cos(y))*cos(y))/k^8 - (2880*exp(k*cos(y))*(7*cos(y)^2 - 1))/k^7 - (exp(k*cos(y))*(cos(y)^2 - 1)^4)/k - (48*exp(k*cos(y))*(35*cos(y)^4 - 30*cos(y)^2 + 3))/k^5 + (960*exp(k*cos(y))*cos(y)*(7*cos(y)^2 - 3))/k^6 + (8*exp(k*cos(y))*cos(y)*(cos(y)^2 - 1)^3)/k^2 - (8*exp(k*cos(y))*(cos(y)^2 - 1)^2*(7*cos(y)^2 - 1))/k^3 + (48*exp(k*cos(y))*cos(y)*(7*cos(y)^4 - 10*cos(y)^2 + 3))/k^4;
                otherwise
                    warning('falling back to numerical integration')
                    cdf = @(y) integral(pdf, 0, y);
            end
            
            u = rand(1,n)*cdf(pi); % consider the fact that pdf is not normalized
            
            lower = zeros(1,n);
            upper = pi * ones(1,n);
            for j=1:52
                middle = (lower+upper)/2;
                cdfMiddle = cdf(middle);
                upper(cdfMiddle>u) = middle(cdfMiddle>u);
                lower(cdfMiddle<=u) = middle(cdfMiddle<=u);
            end
            phi = middle;
            W = cos(phi);
            for i=1:n
                % Step 3 from Wood
                V = mvnrnd(zeros(this.dim-1,1), eye(this.dim-1,this.dim-1));
                V = V'/norm(V);
                % we differ from Wood's algorithm here and sample a
                % distribution with mu = [1,0,...,0] rather than mu =
                % [0,...,0,1]
                X = [W(i); sqrt(1-W(i)^2)*V]; 
                
                s(:,i) = Q*X;
            end
        end
        
        function e = entropy(this)
            % Computes the entropy of the distribution using an analytical
            % solution.
            %
            % Returns:
            %   e (scalar)
            %       the computed entropy
            arguments
                this (1,1) VMFDistribution
            end
            p = this.dim;
            % https://www.wolframalpha.com/input/?i=diff%281%2F%28kappa%5E%28p%2F2-1%29%2F%282*pi%29%5E%28p%2F2%29%2Fbesseli%28p%2F2-1%2C+kappa%29%29%2C+kappa%29
            Cinvdiff = 2^(p/2-1) *pi^(p/2)* this.kappa^(1-p/2) * (besseli(p/2-2,this.kappa)+besseli(p/2,this.kappa))+(1-p/2)*(2*pi)^(p/2)* this.kappa^(-p/2) *besseli(p/2-1,this.kappa);
            e = -log(this.C) - this.kappa * this.C * Cinvdiff;
        end
        
        function h = hellingerDistance(this, other)
            % Analytically calculates the Hellinger distance to another
            % VMF distribution.
            %
            % Parameters:
            %   other (VMFDistribution)
            %       distribution to compare to
            % Returns:
            %   dist (scalar)
            %       hellinger distance of this distribution to other distribution            
            assert(isa(other,'VMFDistribution'));
            assert(this.dim==other.dim);
            d = this.dim;
            kx = this.kappa;
            ky = other.kappa;
            ks = norm(kx/2*this.mu + ky/2*other.mu);
            h = 1- sqrt(kx^(d/2-1)*ky^(d/2-1))/ks^(d/2-1) * besseli(d/2-1,ks)/sqrt(besseli(d/2-1,kx)*besseli(d/2-1,ky));
        end
        
        function r = moment(this)
            % Returns the mean resultant vector
            arguments
                this (1,1) VMFDistribution
            end
            r = VMFDistribution.Ad(this.dim, this.kappa)*this.mu;
        end

        function distShifted = shift(this, offsets)
            % There is no true shifting for the hypersphere. This is a function for compatibility and only works when mu is [0,0,...,1].
            arguments
                this (1,1) VMFDistribution
                offsets (:,1) double {mustBeNonempty}
            end
            assert(isequal(this.mu,[zeros(this.dim-1,1);1]), 'There is no true shifting for the hypersphere. This is a function for compatibility and only works when mu is [0,0,...,1].');
            distShifted = this;
            distShifted.mu = offsets;
        end
    end
    
    methods (Static)
        function V = fit(samples, weights)
            % Fits VMF parameters to a set of samples.
            % Parameters:
            %   samples (d x n matrix)
            %       matrix that contains one sample per column
            %   weights (1 x n row vector)
            %       weight for each sample            
            % Returns:
            %   V (VMFDistribution)
            %       the MLE estimate for a VMF distribution given the
            %       samples
            %
            % Sra, S. 
            % A Short Note on Parameter Approximation for von Mises--Fisher
            % Distributions: And a Fast Implementation of Is (x)
            % Computational Statistics, Springer-Verlag, 2012, 27, 177-190

            d = size(samples,1);    
            n = size(samples,2);
            if nargin<2
                weights = ones(1,n)/n;
            end
            assert(size(weights,1)==1, 'weights needs to be a row vector');
            assert(size(weights,2)==n, 'number of weights and samples needs to match');
            assert(abs(sum(weights)-1) < 0.001, 'weights must sum to 1'); %check normalization
            
            vecSum = sum(samples.*repmat(weights,d,1),2);
            mu_ = vecSum/norm(vecSum);
            Rbar = norm(vecSum);
            kappa_= VMFDistribution.AdInverse(d, Rbar);
            
            V = VMFDistribution(mu_, kappa_);
        end
        
        function V = fitScoreBased(samples, weights)
            % Fits VMF parameters to a set of samples.
            % Parameters:
            %   samples (d x n matrix)
            %       matrix that contains one sample per column
            %   weights (1 x n row vector)
            %       weight for each sample            
            % Returns:
            %   V (VMFDistribution)
            %       the MLE estimate for a VMF distribution given the
            %       samples
            %
            % Mardia, K. V.; Kent, J. T. & Laha, A. K. 
            % Score Matching Estimators for Directional Distributions 
            % arXiv preprint arXiv:1604.08470, 2016
            % Sec. 7.1         
            
            dim = size(samples,1);    
            n = size(samples,2);
            if nargin<2
                weights = ones(1,n)/n;
            else
                error('weighted samples not yet supported'); %todo add support for weights
            end
            assert(size(weights,1)==1, 'weights needs to be a row vector');
            assert(size(weights,2)==n, 'number of weights and samples needs to match');
            assert(abs(sum(weights)-1) < 0.001, 'weights must sum to 1'); %check normalization
            
            vecSum = sum(samples.*repmat(weights,dim,1),2);
            mu_ = vecSum/norm(vecSum);
            
            % use score based method to get kappa
            tmp = VMFDistribution(mu_, 1);
            R = tmp.getRotationMatrix();
            Y = samples' * R;
            
            W = 1 - 1/n * sum(Y(:,1).^2);
            d = (dim-1)/n*sum(Y(:,1)); %/n is missing in Mardia's paper
            kappa_ =d/W;
            
            V = VMFDistribution(mu_, kappa_);            
        end
        
        function V = fromMoment(m)
            % Estimates paramters from mean resultant vector
            
            assert(size(m,2) == 1, 'mu must be a column vector');
            assert(size(m,1) >= 2, 'mu must be at least 2 for the circular case');
            
            mu_ = m/norm(m);
            Rbar = norm(m);
            kappa_= VMFDistribution.AdInverse(length(m), Rbar);
            
            V = VMFDistribution(mu_, kappa_);
        end
        
        function result = Ad(d, kappa)
            % Ratio of bessel functions.
            % Parameters:
            %   d (scalar)
            %       dimension
            %   kappa (scalar)
            %       concentration
            % Returns:
            %   result (scalar)
            %       computed ratio
            %
            % Attention: Ad is defined differently from
            % besselratio/besselratioinverse
            result = besseli(d/2, kappa, 1)/besseli(d/2-1, kappa, 1);
        end
        
        function result = AdInverse(d, x)
            % Computed the inverse of Ad.
            % Parameters:
            %   d (scalar)
            %       dimension
            %   x (scalar)
            %       value of Ad(d, kappa)
            % Returns:
            %   result (scalar)
            %       computed inverse
            %
            % Sra, S. 
            % A Short Note on Parameter Approximation for von Mises--Fisher 
            % Distributions: And a Fast Implementation of Is(x)
            % Computational Statistics, Springer-Verlag, 2012, 27, 177-190
            
            % approximation by Banaerjee et al.
            kappa_ = x*(d-x^2)/(1-x^2);

            % newton approximation by Sra 2012
            newtonStep = @(kappa) kappa - (VMFDistribution.Ad(d, kappa)-x)/(1- VMFDistribution.Ad(d, kappa)^2 - (d-1)/kappa*VMFDistribution.Ad(d, kappa));
            
            maxSteps = 20;
            epsilon = 1E-7; % stopping condition
            
            for i=1:maxSteps
                kappaOld = kappa_;
                kappa_ = newtonStep(kappa_);
                if abs(kappa_ - kappaOld) < epsilon
                    break;
                end
            end
            
            result = kappa_;
        end
    end
end