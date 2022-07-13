classdef LinBoundedDiracDistribution < AbstractLinBoundedDistribution & AbstractDiracDistribution
    methods
        function this = LinBoundedDiracDistribution(boundD_, d_, w_)
            % Constructor, w_ is optional
            %
            % Parameters:
            %   boundD (scalar)
            %       number of circular dimensions
            %   d_ ((boundD_+linD_) x L)
            %       Dirac locations - first bounded then linear dimensions
            %   w_ (1 x L)
            %       weights for each Dirac
            arguments
                % The order is like that to make w_ optional
                boundD_ (1,1) double {mustBeInteger,mustBePositive}
                d_ (:,:) double
                w_ (1,:) double = ones(1,size(d_,2))/size(d_,2); % All Diracs have equal weights by default
            end
            linD_ = size(d_,1) - boundD_; 
            this@AbstractDiracDistribution(d_, w_);
            this.linD = linD_;
            this.boundD = boundD_;
        end
        
        function dist = marginalizePeriodic(this)
            % Marginalizes the periodic dimensions, returns a
            % LinearDiracDistribution of the remaining linear dimensions
            %
            % Returns:
            %   dist (LinearDiracDistribution)
            %       marginal distribution
            arguments
                this (1,1) LinBoundedDiracDistribution
            end
            dist = LinearDiracDistribution(this.d(this.boundD+1:end,:), this.w);
        end
        function m = mean(~) %#ok<STOUT> 
            error('Cannot provide mean without knowledge about the topology of the bounded domain.')
        end
    end
    methods (Sealed)
        function m = linearMean(this)
            arguments
                this (1,1) LinBoundedDiracDistribution
            end
            m = this.marginalizePeriodic().mean();
        end
        function C = linearCovariance(this)
            % Computes covariance of linear dimensions
            %
            % Returns:
            %   C (linD x linD)
            %       covariance matrix
            arguments
                this (1,1) LinBoundedDiracDistribution
            end
            C = this.marginalizePeriodic().covariance();
        end 
    end
end