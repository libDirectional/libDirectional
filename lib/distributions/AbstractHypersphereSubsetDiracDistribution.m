classdef (Abstract) AbstractHypersphereSubsetDiracDistribution < AbstractDiracDistribution & AbstractHypersphereSubsetDistribution
    % Dirac distribution on the hypersphere with dirac positions d and
    % weights w.
    
    methods
        function this = AbstractHypersphereSubsetDiracDistribution(d_, w_)
            arguments
                d_ (:,:) double {mustBeNonempty}
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2);
            end
            assert ( max(abs(vecnorm(d_)-1)) < 1E-5,'DiracLocation:notOnSphere');
            this@AbstractDiracDistribution(d_, w_);
        end

        function C = moment(this)
            C = this.d.*this.w*this.d';
        end
        
        function result = entropy(this)
            % Calculates the entropy analytically 
            %
            % Returns:
            %   result (scalar)
            %       entropy of the distribution
            warning('ENTROPY:DISCRETE','entropy is not defined in a continous sense')
            % return discrete entropy
            % The entropy only depends on the weights!
            result = -sum(this.w.*log(this.w));
        end

        function result = integral(this)
            arguments
                this (1,1) AbstractHypersphereSubsetDiracDistribution
            end
            result = integral@AbstractDiracDistribution(this);
        end

        function s = sample(this, n)
            s = sample@AbstractDiracDistribution(this, n);
        end
    end
end