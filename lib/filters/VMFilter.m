classdef VMFilter < AbstractCircularFilter
    % A filter based on the VM distribution.
    %
    % based on
    % M. Azmani, S. Reboul, J.-B. Choquel, and M. Benjelloun, "A recursive
    % fusion filter for angular data" in 2009 IEEE International Conference
    % on Robotics and Biomimetics (ROBIO), Dec. 2009, pp. 882-887.
    %
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Recursive Bayesian Filtering in Circular State Spaces
    % arXiv preprint: Systems and Control (cs.SY), January 2015.    
    %
    % Igor Gilitschenski, Gerhard Kurz, Uwe D. Hanebeck,
    % Non-Identity Measurement Models for Orientation Estimation Based on Directional Statistics
    % Proceedings of the 18th International Conference on Information Fusion (Fusion 2015),
    % Washington D. C., USA, July 2015.
    
    properties
        vm VMDistribution
    end
    
    methods
        
        function this = VMFilter()
            % Constructor   
            vm_ = VMDistribution(0,1);            
            this.setState(vm_);
        end
        
        function setState(this, vm_)
            % Sets the current system state
            %
            % Parameters:
            %   vm_ (VMDistribution)
            %       new state
            assert(isa(vm_, 'VMDistribution'));
            this.vm = vm_;
        end
        
        function predictIdentity(this, vmSys)
            % Predicts assuming identity system model, i.e.,
            % x(k+1) = x(k) + w(k)    mod 2pi,
            % where w(k) is additive noise given by vmSys.
            %
            % Parameters:
            %   vmSys (VMDistribution)
            %       distribution of additive noise
            assert(isa(vmSys, 'VMDistribution'));
            this.vm = this.vm.convolve(vmSys);
        end
        
        function predictNonlinear(this, f, vmSys) 
            % Predicts assuming a nonlinear system model, i.e.,
            % x(k+1) = f(x(k)) + w(k)    mod 2pi,
            % where w(k) is additive noise given by vmSys.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2pi) to [0,2pi)
            %   vmSys (VMDistribution)
            %       distribution of additive noise
            assert(isa(vmSys, 'VMDistribution'));
            assert(isa(f,'function_handle'));
            
            wd = this.vm.toDirac5();
            wdF = wd.applyFunction(f);
            vmF = wdF.toVM();
            this.vm = vmF.convolve(vmSys);
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
            
            wd = this.vm.toDirac5();
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
            this.vm = wdPosterior.toVM();
        end        
        
        function updateIdentity(this, vmMeas, z)
            % Updates assuming identity measurement model, i.e.,
            % z(k) = x(k) + v(k)    mod 2pi,
            % where v(k) is additive noise given by vmMeas.
            %
            % Parameters:
            %   vmMeas (VMDistribution)
            %       distribution of additive noise
            %   z (scalar)
            %       measurement in [0, 2pi)
            assert(isa(vmMeas, 'VMDistribution'));
            
            muWnew = mod(z - vmMeas.mu,2*pi);
            vmMeasShifted = VMDistribution(muWnew, vmMeas.kappa);
            this.vm = this.vm.multiply(vmMeasShifted);
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
            
            wd = this.vm.toDirac5();
            wdNew = wd.reweigh(@(x) likelihood(z,x));
            this.vm = wdNew.toVM();
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
            samples = this.vm.sample(n);
            wd = WDDistribution(samples);
            wdNew = wd.reweigh(@(x) likelihood(z,x));
            this.vm = wdNew.toVM();
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
            if nargin<4
                tau = 0.02; 
            end
            lambda = 1;
            steps = 0;
            while lambda> 0
                wd = this.vm.toDirac3();
                l = (arrayfun(@(x) likelihood(z,x), wd.d));
                lmin = min(l);
                lmax = max(l);
                if lmax == 0
                    warning('progressive update failed because likelihood is 0 everwhere')
                    return
                end
                wmin = min(wd.w);
                assert(wmin>0);
                wmax = max(wd.w);
                currentLambda = min(log(tau*wmax/wmin)/log(lmin/lmax), lambda);
                if currentLambda <= 0
                    warning('progressive update with given threshold impossible')
                    currentLambda = 0.001;
                end
                wdNew = wd.reweigh(@(x) likelihood(z,x).^currentLambda);
                this.vm = wdNew.toVM();
                lambda = lambda - currentLambda;
                steps = steps + 1;
            end
        end
        
        function updateNonlinearProgressiveJannik(this, likelihood, z, tau)
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
            if nargin<4
                tau = 0.02; 
            end
            lambda = 1;
            steps = 0;
            vmTemp = VMDistribution(this.vm.mu, 0);
            while lambda> 0
                wd = vmTemp.toDirac3();
                l = (arrayfun(@(x) likelihood(z,x).*this.vm.pdf(x), wd.d));
                lmin = min(l);
                lmax = max(l);
                if lmax == 0
                    warning('progressive update failed because likelihood is 0 everwhere')
                    return
                end
                wmin = min(wd.w);
                assert(wmin>0);
                wmax = max(wd.w);
                currentLambda = min(log(tau*wmax/wmin)/log(lmin/lmax), lambda);
                if currentLambda <= 0
                    warning('progressive update with given threshold impossible')
                    currentLambda = 0.001;
                end
                wdNew = wd.reweigh(@(x) likelihood(z,x).^currentLambda.*this.vm.pdf(x).^currentLambda);
                vmTemp = wdNew.toVM();
                lambda = lambda - currentLambda;
                steps = steps + 1;
            end
            this.vm = vmTemp;
        end
        
        function updateNonlinearShift(this, z, h, vmMeas)
            % Update step approximating h(x) by a shift.
            %
            %  Description:
            %    The measurement model assumed by this update step is given
            %    by a shift approximation of h. That is, the true model
            %      z = h(x) + v
            %    is approximated by the simpler model
            %      z = x + s + v .
            %
            %  Parameters:
            %    z - Observation.
            %    h - Observation Function.
            %    vmMeas - Distribution of observation noise.
            
            assert(isa(vmMeas, 'VMDistribution'));
            assert(isa(h,'function_handle'));
            
            % Compute approximated shift
            s = mod(h(this.vm.mu)-this.vm.mu,2*pi);
            
            % Incorporate Measurement
            muNew = mod(z-(s+vmMeas.mu),2*pi);
            
            % Perform measurement update.
            vmObsShifted = VMDistribution(muNew, vmMeas.kappa);
            this.vm = this.vm.multiply(vmObsShifted);
        end
        
        function updateNonlinearStatisticalShift(this, z, h, vmMeas, samplingMethod, samplingParameter)
            % Update step approximating h(x) by a shift.
            %
            %  Description:
            %    The measurement model assumed by this update step is given
            %    by a shift approximation of h. That is, the true model
            %      z = h(x) + v
            %    is approximated by the simpler model
            %      z = x + s + v .
            %
            %  Parameters:
            %    z - Observation.
            %    h - Observation Function.
            %    vmMeas - Distribution of observation noise.
            %    samplingMethod - May be 'dirac3' (default) or 'dirac5'
            %    samplingParameter - Parameter for dirac5 sampling (0.5 is
            %    default).
            
            assert(isa(vmMeas, 'VMDistribution'));
            assert(isa(h,'function_handle'));
            
            % Deterministically Sample system State
            if nargin < 5
                wdState = this.vm.toDirac3;
            elseif strcmp(samplingMethod,'dirac5')
                if nargin < 6
                    samplingParameter = 0.5;
                end
                wdState = this.vm.toDirac5(samplingParameter);
            else % Default
                wdState = this.vm.toDirac3;
            end
            
            % Compute sample based approximation of shift
            s = atan2(sin(h(wdState.d)-wdState.d)*wdState.w', ...
                cos(h(wdState.d)-wdState.d)*wdState.w');
            
            %s = mod(h(this.vm.mu)-this.vm.mu,2*pi);
            
            % Incorporate Measurement
            muNew = mod(z-(s+vmMeas.mu),2*pi);
            
            % Perform measurement update.
            vmObsShifted = VMDistribution(muNew, vmMeas.kappa);
            this.vm = this.vm.multiply(vmObsShifted);
        end
        
        function updateNonlinearCorrectedNoise(this, z, h, vmMeas, samplingMethod, samplingParameter)
            % Update step approximating h(x)+v_1 by x+v_2
            %
            %  Description:
            %    The measurement model assumed by this update step is given
            %    by a shift approximation of h. That is, the true model
            %      z = h(x) + v_1
            %    is approximated by the simpler model
            %      z = x + v_2 ,
            %    where v_2 is also assumed to be a von Mises distributed RV.
            %    The resulting distribution of v_2 will have at most the same
            %    concentration parameter as v_1.
            %
            %  Parameters:
            %    z - Observation.
            %    h - Observation Function.
            %    vmMeas - Distribution of observation noise.
            %    samplingMethod - May be 'dirac3' (default) or 'dirac5'
            %    samplingParameter - Parameter for dirac5 sampling (0.5 is
            %    default).
            
            assert(isa(vmMeas, 'VMDistribution'));
            assert(isa(h,'function_handle'));           
            
            % Deterministically Sample system State
            if nargin < 5
                samplingMethod = 'dirac3';
            end
            
            % Choose sampling method.
            if strcmp(samplingMethod,'dirac5')
                if nargin < 6
                    samplingParameter = 0.5;
                end
                wdState = this.vm.toDirac5(samplingParameter);
                wdNoise = vmMeas.toDirac5(samplingParameter);
            elseif strcmp(samplingMethod,'dirac3') % Default 'dirac3'
                wdState = this.vm.toDirac3;
                wdNoise = vmMeas.toDirac3;
            else
                error('Unknown sampling method: %s', samplingMethod);
            end
            
            % Measurement difference function
            hDiff = @(x,v) h(x) - x + v;
            
            k = 1;
            newDiracs = zeros(1,numel(wdState.d)*numel(wdNoise.d));
            newWeights= zeros(1,numel(wdState.d)*numel(wdNoise.d));
            for i=1:numel(wdState.d)
                for j=1:numel(wdNoise.d)
                    newDiracs(k) = hDiff(wdState.d(i),wdNoise.d(j));
                    newWeights(k) = wdState.w(i)*wdNoise.w(j);
                    k=k+1;
                end
            end
            wdNewNoise = WDDistribution(newDiracs,newWeights);
            vmNewNoise = wdNewNoise.toVM;
            
            % Now perform typical identity based update.
            muWnew = mod(z - vmNewNoise.mu,2*pi);
            vmMeasShifted = VMDistribution(muWnew, vmNewNoise.kappa);
            this.vm = this.vm.multiply(vmMeasShifted);
        end
        
        function vm = getEstimate(this)
            % Return current estimate 
            %
            % Returns:
            %   vm (VMDistribution)
            %       current estimate
            vm = this.vm;
        end
        
        function likelihoodVal=associationLikelihood(this,likelihood)
            % see Florian Pfaff, Kailai Li, and Uwe D. Hanebeck,
            % Association Likelihoods for Directional Estimation
            % Proceedings of the 2019 IEEE International Conference on
            % Multisensor Fusion and Integration for Intelligent Systems (MFI 2019),
            % Taipei, Republic of China, May, 2019.
            vmEst = this.getEstimate.multiply(likelihood);
            likelihoodVal = besseli(0, vmEst.kappa) ...
                / (2 * pi * besseli(0, this.getEstimate.kappa)...
                * besseli(0, likelihood.kappa));
        end
    end
    
end

