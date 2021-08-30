classdef WNFilter < AbstractCircularFilter
    % A filter based on the WN distribution.
    %   
    % based on
    % G. Kurz, I. Gilitschenski, and U. D. Hanebeck, “Recursive nonlinear
    % filtering for angular data based on circular distributions,” in Proceedings
    % of the 2013 American Control Conference (ACC 2013),
    % Washington D.C., USA, Jun. 2013.
    %
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Nonlinear Measurement Update for Estimation of Angular Systems Based on Circular Distributions
    % Proceedings of the 2014 American Control Conference (ACC 2014), 
    % Portland, Oregon, USA, June 2014.
    %
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Recursive Bayesian Filtering in Circular State Spaces
    % arXiv preprint: Systems and Control (cs.SY), January 2015.
        
    properties
        wn (1,1) WNDistribution = WNDistribution(0,1) % Initial Estimate 
    end
    
    methods
        
        function setState(this, wn_)
            % Sets the current system state
            %
            % Parameters:
            %   wn_ (WNDistribution)
            %       new state
            assert(isa(wn_, 'WNDistribution'));
            this.wn = wn_;
        end
        
        function predictIdentity(this, wnSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by wnSys.
            %
            % Parameters:
            %   wnSys (WNDistribution)
            %       distribution of additive noise
            assert(isa(wnSys, 'WNDistribution'));
            this.wn = this.wn.convolve(wnSys);
        end
        
        function predictNonlinear(this, f, wnSys) 
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by wnSys.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi)
            %   wnSys (WNDistribution)
            %       distribution of additive noise
            assert(isa(wnSys, 'WNDistribution'));
            assert(isa(f,'function_handle'));
            
            wd = this.wn.toDirac5();
            wdF = wd.applyFunction(f);
            wnF = wdF.toWN();
            this.wn = wnF.convolve(wnSys);
        end
        
        function predictNonlinearNonAdditive(this, f, samples, weights)
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k), w(k))    mod 2pi,
            % where w(k) is non-additive noise given by samples and weights.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) x W to [0,2pi) (W is the space
            %       containing the noise samples)
            %   samples (d x n matrix)
            %       n samples of the noise as d-dimensional vectors
            %   weights (1 x n vector)
            %       weight of each sample
           
            % (samples, weights) are discrete approximation of noise
            assert(size(weights,1) == 1, 'weights most be row vector')
            assert(size(samples,2) == size(weights,2), 'samples and weights must match in size');
            assert(isa(f,'function_handle'));
            
            weights = weights/sum(weights); % ensure normalization of weights
            
            wd = this.wn.toDirac5();
            L = length(wd.d);
            Lw = length(samples);
            d = zeros(1,L*Lw);
            w = zeros(1,L*Lw);
            for j=1:L
                for l=1:Lw
                    w(1,j+L*(l-1)) = wd.w(j) * weights(l);
                    d(1,j+L*(l-1)) = f(wd.d(j), samples(:,l));
                end
            end
            wdPosterior = WDDistribution(d,w);
            this.wn = wdPosterior.toWN();
        end

        function updateIdentity(this, wnMeas, z) 
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by wnMeas.
            %
            % Parameters:
            %   wnMeas (WNDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            assert(isa(wnMeas, 'WNDistribution'));
            
            muWnew = mod(z - wnMeas.mu,2*pi);
            wnMeasShifted = WNDistribution(muWnew, wnMeas.sigma);
            this.wn = this.wn.multiplyVM(wnMeasShifted);
        end
        
        function updateNonlinear(this, likelihood, z) 
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            %
            % This method uses a simple reweighting algorithm.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi) to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            assert(isa(likelihood,'function_handle'));
            
            wd = this.wn.toDirac5();
            wdNew = wd.reweigh(@(x) likelihood(z,x));
            this.wn = wdNew.toWN();
        end
        
        function updateNonlinearParticle(this, likelihood, z) 
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            %
            % This method uses a nondeterministic Monte Carlo algorithm.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi) to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            assert(isa(likelihood,'function_handle'));
            
            n = 100;
            samples = this.wn.sample(n);
            wd = WDDistribution(samples);
            wdNew = wd.reweigh(@(x) likelihood(z,x));
            this.wn = wdNew.toWN();
        end
        
        function updateNonlinearProgressive(this, likelihood, z, tau) 
            % Updates assuming nonlinear measurement model given by a
            % likelihood function likelihood(z,x) = f(z|x), where z is the
            % measurement. The function can be created using the
            % LikelihoodFactory.
            %
            % This method uses a deterministic progressive algorithm.
            % 
            % Parameters:
            %   likelihood (function handle)
            %       function from Z x [0,2pi) to [0, infinity), where Z is
            %       the measurement space containing z
            %   z (arbitrary)
            %       measurement
            %   tau (scalar between 0 and 1)
            %       parameter controlling the progression step size
            assert(isa(likelihood,'function_handle'));
            
            if nargin<4
                tau = 0.02; 
            end
            lambda = 1;
            steps = 0;
            while lambda> 0
                wd = this.wn.toDirac5();
                l = (arrayfun(@(x) likelihood(z,x), wd.d));
                lmin = min(l);
                lmax = max(l);
                if lmax == 0
                    warning('progressive update failed because likelihood is 0 everwhere')
                    return
                end
                wmin = min(wd.w);
                assert(wmin > 0);
                wmax = max(wd.w);
                currentLambda = min(log(tau*wmax/wmin)/log(lmin/lmax), lambda);
                if currentLambda <= 0
                    warning('progressive update with given threshold impossible')
                    currentLambda = 0.001;
                end
                wdNew = wd.reweigh(@(x) likelihood(z,x).^currentLambda);
                this.wn = wdNew.toWN();
                lambda = lambda - currentLambda;
                steps = steps + 1;
            end
        end
        
        function wn = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   wn (WNDistribution)
            %       current estimate
            wn = this.wn;
        end
        
    end
    
end

