classdef (Abstract) AbstractDummyFilter < AbstractFilter
    properties (SetAccess = protected)
        dist 
    end
    methods
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
            hfd=this.dist;
        end    
        
        function est = getPointEstimate(this)
            est = this.dist.sample(1);
        end        
    end
end