classdef (Abstract) AbstractGridDistribution < AbstractDistribution
    
    properties
        gridValues double {mustBeNonnegative}
        enforcePdfNonnegative logical = true
    end
    methods (Abstract)
        getManifoldSize(this);
    end
    methods
        function pdf(~, ~)
            % If no interpolation is defined in the subclass, return error
            error('PDF:UNDEFINED', 'Pdf not defined');
        end
        
        function gd = multiply(this, other)
            arguments
                this AbstractGridDistribution
                other AbstractGridDistribution
            end
            assert(this.enforcePdfNonnegative == other.enforcePdfNonnegative);
            gd = this;
            gd.gridValues = gd.gridValues.*other.gridValues;
            warnStruct = warning('off', 'Normalization:notNormalized');
            gd = gd.normalize;
            warning(warnStruct);
        end
        
        function gd = normalize(this, tol)
            arguments
                this AbstractGridDistribution
                tol (1,1) double = 1e-4
            end
            integral = this.getManifoldSize * mean(this.gridValues);
            gd = this; % Asssign here so we can use return in the else branch
            if any(this.gridValues < 0)
                warning('Normalization:negative', 'There are negative values. This usualy points to a user error.');
            elseif abs(integral) < 1e-200 % Tolerance has to be that low to avoid unnecessary errors on multiply
                error('Normalization:almostZero', 'Sum of grid vallues is too close to zero, this usually points to a user error.');
            elseif abs(integral-1) > tol
                warning('Normalization:notNormalized', 'Grid values apparently do not belong to normalized density. Normalizing...');
            else % Density is normalized. Do not change anything
                return
            end
            gd.gridValues = gd.gridValues/integral;
        end

                
    end
end

