classdef (Abstract) AbstractDiracDistribution < AbstractDistribution
    properties
        d (:,:) double
        w (1,:) double
    end
    
    methods
        function this = AbstractDiracDistribution(d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   d_ (dim x L)
            %       Dirac locations
            %   w_ (1 x L)
            %       weights for each Dirac
            arguments
                d_ (:,:) double {mustBeNonempty}
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            this.dim = size(d_,1);
            if size(d_,1)>size(d_,2)
                warning('Not even one one Dirac per dimension. If this warning is unxpected, verify d_ is shaped correctly.');
            end
            assert(size(d_,2) == size(w_,2),'Number of Dircas and weights must match.');
            this.d = d_;
            if ~(abs(sum(w_)-1)<1e-10)
                warning('Dirac:WeightsUnnormalized','Sum of the weights is not 1. Normalizing to 1')
                this.w = w_/sum(w_);
            else
                this.w = w_;
            end
        end
                
        function dist = applyFunction(this,f)
            % Apply a function f(x) to each Dirac component and obtain its new position
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi)^dim to [0,2*pi)^dim
            % Returns:
            %   twd (ToroidalWDDistribution)
            %       distribution with new Dirac locations (and same
            %       weights as before)
            arguments
                this (1,1) AbstractDiracDistribution
                f (1,1) function_handle
            end
            d_ = zeros(size(this.d));
            for i=1:size(this.d,2)
                d_(:,i) = f(this.d(:,i));
            end
            dist = this;
            dist.d = d_;
        end
        
        function dist = reweigh(this, f)
            % Uses a function f(x) to calculate the weight of each Dirac
            % component. The new weight is given by the product of the old 
            % weight and the weight obtained with f. Restores normalization
            % afterwards.
            %
            % Parameters:
            %   f (function handle)
            %       function from [0,2*pi)^dim to [0, infinity)
            %       (needs to support vectorized inputs, i.e., dim x n matrices)
            % Returns:
            %   wd (some DiracDistribution)
            %       distribution with new weights and same Dirac locations
            arguments
                this (1,1) AbstractDiracDistribution
                f (1,1) function_handle
            end
            wNew = f(this.d);
            assert(isequal(size(wNew),[1,size(this.d,2)]),'Function returned wrong number of outputs.');
            assert(all(wNew >= 0));
            assert(sum(wNew) > 0);
            
            dist = this;
            dist.w = wNew.*this.w;
            dist.w = dist.w/sum(dist.w);
        end

        function s = sample(this, n)
            % Obtain n samples from the distribution
            %
            % Parameters:
            %   n (scalar)
            %       number of samples
            % Returns:
            %   s (dim x n)
            %       one sample per column
            arguments
                this (1,1) AbstractDiracDistribution
                n (1,1) double {mustBeInteger,mustBePositive}
            end
            ids = discretesample(this.w,n);
            s = this.d(:,ids);
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            arguments
                this (1,1) AbstractDiracDistribution
            end
            warning('ENTROPY:DISCRETE','entropy is not defined in a continous sense')
            % return discrete entropy
            % The entropy only depends on the weights!
            result = -sum(this.w.*log(this.w));
        end

        function result = integral(this)
            % Calculate integral to check normalization
            % (should always be 1)
            % Returns:
            %   result (scalar)
            %       integral over hypersphere surface of pdf (uses
            %       approximation, not very accurate for higher dimensions)
            %  Cannot be sealed because additional arguments
            % (integration limits) for 1-D case !!!
            arguments
                this (1,1) AbstractDiracDistribution
            end
            result = sum(this.w);
        end
    end
    
    methods (Sealed)
        function logLikelihood(~, ~)
            error('PDF:UNDEFINED', 'not supported');
        end     

        function pdf(~, ~)
            % Placeholder, pdf does not exist for wrapped Dirac distributions
            error('PDF:UNDEFINED', 'pdf is not defined')
        end    
        
        function integralNumerical(~, ~)
            error('PDF:UNDEFINED', 'not supported')
        end        
        
        function trigonometricMomentNumerical(~,~)
            % Disable numerical calculation of angular moments since it relies on the pdf
            error('PDF:UNDEFINED', 'not supported');
        end   
        
        function sampleMetropolisHastings(~, ~)
            % Disable sampling algorithm relying on pdf
            error('PDF:UNDEFINED', 'not supported');
        end        

        function squaredDistanceNumerical(~, ~)
            error('PDF:UNDEFINED', 'not supported');
        end
        
        function kldNumerical(~, ~)
            error('PDF:UNDEFINED', 'not supported');
        end       

        function m = mode(this, relTol)
            arguments
                this (1,1) AbstractDiracDistribution
                relTol (1,1) double = 0.001
            end
            % sample with maximum weight
            [highestVal,ind] = max(this.w);
            if (highestVal/numel(this.w))<(1+relTol)
                warning('The samples may be equally weighted, .mode is likely to return a bad result.')
            end
            m = this.d(:,ind);
        end

        function modeNumerical(~)
            error('PDF:UNDEFINED', 'not supported');
        end

        function entropyNumerical(~)
            error('PDF:UNDEFINED', 'not supported');
        end
    end
end

