classdef (Abstract) AbstractDummyFilter < AbstractFilter
    methods
        function this = AbstractDummyFilter(dim)
            arguments
                dim (1,1) double {mustBeGreaterThanOrEqual(dim,1)}
            end
            this.dim=dim;
            % Do nothing
        end
        
        function setState(this, dist)
            assert(dist.dim==this.dim);
            % Do nothing
        end
        
        function predictIdentity(~, ~)
            % Do nothing
        end
        
        function predictNonlinear(~, ~, ~,~)
            % Do nothing
        end
        
        function predictNonlinearViaTransitionDensity(~, ~,~)
            % Do nothing
        end
        
        function updateIdentity(~, ~, ~)
            % Do nothing
        end

        function updateNonlinear(~, ~, ~)
            % Do nothing
        end
        
        function hfd = getEstimate(this)
            hfd=HypersphericalUniformDistribution(this.dim);
        end
        
        function mean=getEstimateMean(this)
            mean=HypersphericalUniformDistribution(this.dim).sample(1);
        end
        
    end
end