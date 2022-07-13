classdef HypercylindricalDiracDistribution < AbstractHypercylindricalDistribution & LinBoundedDiracDistribution
    % Partially wrapped dirac distribution with Dirac positions d and
    % weights w
    %
    % Gerhard Kurz
    % Directional Estimation for Robotic Beating Heart Surgery
    % Karlsruhe Institute of Technology, 
    % Intelligent Sensor-Actuator-Systems Laboratory, 2015               
    %
    
    methods
        function this = HypercylindricalDiracDistribution(boundD_, d_, w_)
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
            d_(1:boundD_,:) = mod(d_(1:boundD_,:),2*pi);
            this@LinBoundedDiracDistribution(boundD_, d_, w_);
        end
        
        function mu = hybridMoment(this)
            % Calculates mean of [cos(x1), sin(x1), cos(x2), .., cos(x_boundD), sin(x_boundD),
            % x_(boundD+1, x_(boundD+2), ..., x_(linD+boundD)]
            %
            % Returns:
            %   mu (linD + 2*boundD)
            %       Calculates mean of [cos(x1), sin(x1), cos(x2), .., cos(x_boundD), sin(x_boundD),
            %       x_(boundD+1, x_(boundD+2), ..., x_(linD+boundD)]
            arguments
                this (1,1) HypercylindricalDiracDistribution
            end
            S = NaN(this.boundD*2+this.linD,size(this.d,2));
            S(2*this.boundD+1:end, :) = this.d(this.boundD+1:end,:);
            for i=1:this.boundD
                S(2*i-1,:) = cos(this.d(i,:));
                S(2*i,:) = sin(this.d(i,:));
            end
            mu = sum(repmat(this.w,this.linD+2*this.boundD,1).*S,2);
        end
        
        function dist = marginalizeLinear(this)
            % Marginalizes the linear dimensions, returns a
            % HypertoroidalWDDistribution of the remaining circular
            % dimensions
            %
            % Returns:
            %   dist (HyperToroidalWDDistribution)
            %       marginal distribution
            arguments
                this (1,1) HypercylindricalDiracDistribution
            end
            dist = HypertoroidalWDDistribution(this.d(1:this.boundD,:), this.w);
        end
        
        function m = hybridMean(this)
            % 
            %
            % Returns:
            %   m (linD + boundD x 1 column vector)
            %       the mode
            arguments
                this (1,1) HypercylindricalDiracDistribution
            end
            periodic = this.marginalizeLinear();
            linear = this.marginalizePeriodic();
            m = [periodic.meanDirection(); linear.mean()];            
        end

        function result = integral(this)
            arguments
                this (1,1) HypertoroidalWDDistribution
            end
            result = integral@AbstractDiracDistribution(this);
        end
    end
    
     methods (Static)
        function ddist = fromDistribution(dist, nParticles)
            arguments
                dist (1,1) AbstractHypercylindricalDistribution
                nParticles (1,1) {mustBeInteger,mustBePositive}
            end
            ddist = HypercylindricalDiracDistribution(dist.boundD, dist.sample(nParticles), 1/nParticles * ones(1,nParticles));
        end
     end
    
end