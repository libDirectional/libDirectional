classdef (Abstract) AbstractGridFilter < AbstractFilter
    
    properties
        gd
    end
    methods        
        function setState(this, gd_)
            % Sets the current system state
            %
            % Parameters:
            %   gd_ (AbstractDistribution, preferably a GridDistribution)
            %       new state
            arguments
                this (1,1) AbstractGridFilter
                gd_ (1,1) AbstractDistribution
            end
            if ~(isa(gd_, 'AbstractGridDistribution'))
                warning('setState:nonGrid', 'gd_ is not a GridDistribution. Transforming the distrition with a number of coefficients that is equal to that of the filter.');
                gd_ = this.gd.fromDistribution(gd_, numel(this.gd.gridValues), this.gd.enforcePdfNonnegative);
            elseif ~isequal(size(this.gd.gridValues),size(gd_.gridValues))
                warning('setState:noOfGridValuesDiffer', 'New gris has a different number of grid points.')
            end
            this.gd = gd_;
        end
        
        function est = getEstimate(this)
            % Return current estimate
            %
            % Returns:
            %   gd (GridDistribution)
            %       current estimate
            est = this.gd;
        end
        
        function updateNonlinear(this, likelihood, z) %measurement z, likelihood(z,x)=P(Z|X)
            arguments
                this AbstractGridFilter
                likelihood function_handle
                z (:,1) double
            end
            % Assume likelihood can use implicit expansion (for scalars also possible in older Matlab versions)
            gridValsNew=this.gd.gridValues.*reshape(likelihood(z, this.gd.getGrid()),[],1);
            assert(isequal(size(gridValsNew),size(this.gd.gridValues)));
            this.gd.gridValues = gridValsNew;
            warnStruct = warning('off','Normalization:notNormalized');
            this.gd = this.gd.normalize;
            warning(warnStruct);
        end
    end
    
end
