classdef VMDistribution < AbstractCircularDistribution
    % von Mises distribution
    %
    % see Sreenivasa Rao Jammalamadaka and A. SenGupta, "Topics in Circular 
    % Statistics", 2001, Sec. 2.2.4, page 35ff.
    
    properties
        mu
        kappa
    end
    
    methods
        function this = VMDistribution(mu_, kappa_)
            % Constructor
            assert(isscalar(mu_));
            this.mu = mod(mu_,2*pi);
            assert(kappa_>=0);
            this.kappa = kappa_;
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
            p = exp(this.kappa.*cos(xa-this.mu)) ./ (2 * pi * besseli(0,this.kappa));
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
            if n==1
                m = besselratio(0, this.kappa)*exp(1i * n * this.mu);
            elseif n==2
                m = besselratio(0, this.kappa)*besselratio(1, this.kappa)*exp(1i * n * this.mu);
            else
                if this.kappa > 500
                    warning ('result may be inaccurate');
                end
                m = besseli(abs(n),this.kappa,1)/besseli(0, this.kappa,1)*exp(1i * n * this.mu);
            end
        end
                
        function vm = multiply(this, vm2)
            % Multiply two von Mises distribution (exact)
            %
            % Parameters:
            %   vm2 (VMDistribution)
            %       distribution to multiply with
            % Returns:
            %   vm (VMDistribution)
            %       product of this and vm2
            
            assert(isa(vm2,'VMDistribution'));
            
            C = this.kappa * cos(this.mu) + vm2.kappa * cos(vm2.mu);
            S = this.kappa * sin(this.mu) + vm2.kappa * sin(vm2.mu);
            mu_ = mod(atan2(S,C),2*pi);
            kappa_ = sqrt(C^2+S^2);
            vm = VMDistribution(mu_, kappa_);
        end
        
        function vm = convolve(this, vm2)
            % Calculate convolution of two von Mises distributions (approximate)
            %
            % Parameters:
            %   vm2 (VMDistribution)
            %       distribution to convolve with
            % Returns:
            %   vm (VMDistribution)
            %       convolution of this and vm2
            %
            % based on the approximation from 
            % Azmani, M.; Reboul, S.; Choquel, J.-B. & Benjelloun, M. A 
            % Recursive Fusion Filter for Angular Data 
            % IEEE International Conference on Robotics and Biomimetics (ROBIO 2009), 2009, 882-887
            mu_ = mod(this.mu + vm2.mu,2*pi);
            t = besselratio(0, this.kappa)*besselratio(0, vm2.kappa);
            kappa_ = besselratioInverse(0,t,'sraamos');
            vm = VMDistribution(mu_,kappa_);
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            result = -this.kappa * besselratio(0, this.kappa) + log(2*pi*besseli(0,this.kappa));
        end
        
        function samples = sample(this, n)
            % Obtain n samples from the distribution
            % 
            % based on 
            % Devroye, L. Non-Uniform Random Variate Generation Springer, 1986
            % together with the errata 
            % http://luc.devroye.org/errors.pdf
            %
            % originally published in
            % Efficient Simulation of the von Mises Distribution
            % D. J. Best and N. I. Fisher
            % Journal of the Royal Statistical Society. Series C (Applied Statistics)
            % Vol. 28, No. 2 (1979), pp. 152-157
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   samples (1 x n)
            %       n samples on the circle
            samples = zeros(1,n);
            
            r = 1 + sqrt(1+4*this.kappa^2); %tau in Best's paper
            rho = (r - sqrt(2*r))/2/this.kappa;
            s = (1+rho^2)/2/rho; % r in Best's paper
            for j=1:n
                while true
                    U = 2*rand(1)-1; % U in [-1,1]
                    V = rand(1); % V in [0,1]
                    Z = cos(pi*U); % z in Best's paper
                    W = (1 + s*Z)/(s+Z); % f in Best's paper
                    Y = this.kappa*(s-W); % c in Best's paper
                    if log(Y/V)+1-Y>=0 % only checking condition 2 should be sufficient
                        break;
                    end
                end
                samples(1,j) = mod(sign(U)*acos(W) + this.mu, 2*pi); 
            end
        end
    end
    
    methods (Static)
        function vm = fromMoment(m)
            % Obtain a VM distribution from a given first trigonometric moment
            %
            % Parameters:
            %   m (scalar)
            %       first trigonemetric moment (complex number)
            % Returns:
            %   vm (VMDistribution)
            %       VM distribution obtained by moment matching
            mu_ = mod(atan2(imag(m),real(m)),2*pi);
            kappa_ =  besselratioInverse(0,abs(m),'sraamos');
            vm = VMDistribution(mu_,kappa_);
        end
    end
end

