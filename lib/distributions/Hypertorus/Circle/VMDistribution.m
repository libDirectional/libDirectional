classdef VMDistribution < AbstractCircularDistribution
    % von Mises distribution
    %
    % see Sreenivasa Rao Jammalamadaka and A. SenGupta, "Topics in Circular
    % Statistics", 2001, Sec. 2.2.4, page 35ff.
    
    properties
        mu (1,1) double {mustBeReal,mustBeFinite,mustBeNonNan}
        kappa (1,1) double {mustBeNonnegative}
    end
    
    methods
        function this = VMDistribution(mu_, kappa_)
            arguments
                mu_ (1,1) double
                kappa_ (1,1) double
            end
            % Constructor
            this.mu = mod(mu_,2*pi);
            this.kappa = kappa_;
        end
        
        function r = cdf(this, xa, startingPoint)
            % Evaluate cumulative distribution function
            %
            % Parameters:
            %   xa (1 x n)
            %       points where the cdf should be evaluated
            %   startingPoint (scalar)
            %       point where the cdf is zero (starting point can be
            %       [0,2pi) on the circle, default is 0)
            % Returns:
            %   val (1 x n)
            %       cdf evaluated at columns of xa
            arguments
                this (1,1) VMDistribution
                xa (1,:) double
                startingPoint (1,1) double = 0
            end
            
            r = zeros(size(xa));
            for  i=1:size(xa,2)
                r(i) = circVMcdf(xa(i)-this.mu,this.kappa) - circVMcdf(startingPoint-this.mu, this.kappa);
                if r(i)<0
                    r(i) = r(i)+1;
                end
            end
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
            if this.kappa == 0
                m = 0;
                return
            end
            
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
            arguments
                this (1,1) VMDistribution
                vm2 (1,1) VMDistribution
            end            
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
            arguments
                this (1,1) VMDistribution
                vm2 (1,1) VMDistribution
            end
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
            arguments
                this (1,1) VMDistribution
            end
            result = -this.kappa * besselratio(0, this.kappa) + log(2*pi*besseli(0,this.kappa));
        end
        
        function vm = shift(this, angle)
            % Shift distribution by the given angle
            %
            % Parameters:
            %   shiftAngles (scalar)
            %       angle to shift by
            % Returns:
            %   hd (VMDistribution)
            %       shifted distribution
            arguments
                this (1,1) VMDistribution
                angle (1,1) double
            end
            vm = VMDistribution(this.mu+angle, this.kappa); % Mod is done in constructor
        end

        function m = mode(this)
            arguments
                this (1,1) VMDistribution
            end
            m = this.mu;
        end

        function dist = setMode(this, newMode)
            arguments
                this (1,1) VMDistribution
                newMode (1,1) double
            end
            dist = this;
            dist.mu = newMode;
        end
        
        function kld = kld(this, other)
            % Analytically calculates  the Kullback-Leibler divergence to another VM distribution
            % Be aware that the kld is not symmetric.
            %
            % Parameters:
            %   other (AbstractCircularDistribution)
            %       distribution to compare to
            % Returns:
            %   kld (scalar)
            %       kld of this distribution to other distribution
            arguments
                this (1,1) VMDistribution
                other (1,1) VMDistribution
            end
            m1 = this.trigonometricMoment(1);
            kld = log(besseli(0,other.kappa)) - log(besseli(0,this.kappa)) ...
                + this.kappa * real(exp(-1i*this.mu)*m1) ...
                - other.kappa * real(exp(-1i*other.mu)*m1);
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
            arguments
                this (1,1) VMDistribution
                n (1,1) double {mustBePositive, mustBeInteger}
            end
            if this.kappa == 0 % the method by Best does not work is kappa is zero
                cu = CircularUniformDistribution();
                samples = cu.sample(n);
                return
            end
            
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
        
        function [samples, weights] = sampleOptimalQuantization(this, N)
            % Computes optimal quantization of the von Mises distribution.
            %
            % Parameters:
            %   N (scalar)
            %       number of samples
            % Returns:
            %   samples (1 x N)
            %       N samples on the circle
            %   weights (1 x N)
            %       weight for each sample
            %
            % Igor Gilitschenski, Gerhard Kurz, Uwe D. Hanebeck, Roland Siegwart,
            % Optimal Quantization of Circular Distributions
            % Proceedings of the 19th International Conference on Information Fusion (Fusion 2016), Heidelberg, Germany, July 2016.
            arguments
                this (1,1) VMDistribution
                N (1,1) double {mustBePositive, mustBeInteger}
            end
            funWithDeriv = @(x) optimfun(x, this.kappa);
            startval = ((2*pi/(2*N) + (0:N-1)*(2*pi)/N)-pi) * min(1,1/sqrt(this.kappa));
            
            [samples, ~, ~, ~] = ...
                fminunc(funWithDeriv, startval, ...
                optimoptions('fminunc', 'Display','none', 'GradObj', 'on', ...
                'MaxFunEvals', 100000, 'MaxIter', 1000, ...
                'TolFun', 1e-12, 'TolX', 1e-12, 'Algorithm', 'trust-region'));
            
            samples = sort(mod(samples + this.mu,2*pi));
            weights = VMWeights(samples, this.mu, this.kappa);
            
            function [v, g] = optimfun(x,kappa)
                v = VMUtility(x, kappa);
                g = VMUtilityDeriv(x, kappa);
            end
            
            function R = VMUtility( x, kappa )
                % Utility function for computing an optimal L-quantizer
                %
                % Parameters:
                %       x - Positions of the discrete points (must be between -pi and pi).
                %       sigma - Dispersion parameter of Wrapped Normal.
                
                % This is not necessary for the optimization procedure.
                nc = 1/(2*pi*besseli(0,kappa));
                
                VMDensityUnnormalized = @(a) nc*exp(kappa .* cos(a));
                
                R = 0;
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.
                
                % Compute the actual utility.
                for i=1:(numel(x)-2)
                    intfun = @(a) (a-x(i+1)).^2.*VMDensityUnnormalized(a);
                    I = integral(intfun, xBoundary(i),xBoundary(i+1));
                    R = R + I;
                end
            end
            
            function G = VMUtilityDeriv(x, kappa )
                % Gradient of Utility function for an optimal L-quantizer
                %
                % Parameters:
                %       x - Positions of the discrete points (must be between -pi and pi).
                %       kappa - Concentration parameter of von Mises
                
                % This is not necessary for the optimization procedure.
                nc = 1/(2*pi*besseli(0,kappa));

                VMDensity = @(a) nc*exp(kappa .* cos(a));
                
                G = zeros(1,numel(x));
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.
                
                % Compute the actual utility.
                for i=1:(numel(x)-2)
                    intfun = @(a) (-2*a+2*x(i+1)).*VMDensity(a);
                    G(i) = integral(intfun, xBoundary(i),xBoundary(i+1));
                end
            end
            
            function [ W ] = VMWeights(x, mu, kappa )
                nc = 1/(2*pi*besseli(0,kappa));
                VMDensity = @(a) nc*exp(kappa .* cos(a-mu));
                
                W = zeros(1,numel(x));
                
                x = [(x(end)-2*pi) x (x(1)+2*pi)]; % Create x_0 and x_{L+1}.
                xBoundary = (x(1:(end-1)) + x(2:end))/2; % Mean points between the x_i.
                
                for i=1:(numel(x)-2)
                    W(i) = integral(VMDensity, xBoundary(i),xBoundary(i+1));
                end
                
                W = W/sum(W);
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

