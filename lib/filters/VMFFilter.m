classdef VMFFilter < AbstractHypersphericalFilter
    % Filter based von the VMF distribution.
    %
    % Gerhard Kurz, Igor Gilitschenski, Uwe D. Hanebeck,
    % Unscented von Mises-Fisher Filtering
    % IEEE Signal Processing Letters, 2016.    
    
    properties
        state
    end
        
    methods 
        function this = VMFFilter()
            this.setState(VMFDistribution([1;0],1));
        end
        
        function setState(this, state_)
            assert(isa(state_, 'VMFDistribution'));
            this.state = state_;
            this.dim = length(this.state.mu);
        end
        
        function est = getEstimate(this)
           est = this.state; 
        end
        
        function predictIdentity(this, sysNoise)
            assert(isa(sysNoise, 'VMFDistribution'));
            
            %todo how to handle rotation Q
            
            this.state = this.state.convolve(sysNoise);
        end
        
        function predictNonlinear(this, f, sysNoise)
            assert(isa(sysNoise, 'VMFDistribution'));
            assert(isa(f, 'function_handle'));
            
            %todo how to handle rotation Q or non zero-mean sysnoise
            
            stateSamples = this.state.sampleDeterministic();
            
            for i=1:size(stateSamples,2)
                stateSamples(:,i) = f(stateSamples(:,i));
            end
            
            vmf = VMFDistribution.fit(stateSamples);
            this.state = vmf.convolve(sysNoise);
        end
        
        function predictNonlinearArbitraryNoise(this, f, noiseSamples, noiseWeights)
            assert(isa(f, 'function_handle'));
            assert(size(noiseSamples,2) == size(noiseWeights,2));
            assert(size(noiseWeights,1) == 1);
            assert(all(noiseWeights>0));
            
            noiseWeights = noiseWeights/sum(noiseWeights); %normalize weights
            
            stateSamples = this.state.sampleDeterministic();
            k = 1;
            for i=1:size(stateSamples,2)
                for j=1:size(noiseSamples,2)
                    newSamples(:,k) = f(stateSamples(:,i), noiseSamples(:,j));
                    newWeights(k) = noiseWeights(j)/size(stateSamples,2);
                    k = k+1;
                end
            end
            
            this.state = VMFDistribution.fit(newSamples, newWeights);
        end
        
        function updateIdentity(this, z, measNoise)
            assert(isa(measNoise, 'VMFDistribution'));
            assert(size(z,1) == this.dim);
            assert(size(z,2) == 1);
            
            %todo how to handle rotation Q or non zero-mean noise
            
            measNoise.mu = z;
            
            this.state = this.state.multiply(measNoise);
        end
        
        function updateNonlinearProgressive(this)
            %todo later
        end
        
        function mean=getEstimateMean(this)
            mean=this.state.mu;
        end
    end
end
