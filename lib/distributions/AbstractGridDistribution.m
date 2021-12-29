classdef (Abstract) AbstractGridDistribution < AbstractDistribution
    
    properties
        gridValues double {mustBeNonnegative}
        enforcePdfNonnegative logical = true
        gridType char = 'unknown'
    end
    properties (GetAccess = protected, SetAccess = public)
        % Since it may stay empty, prevent people from getting it directly.
        % Use getGrid instead!
        grid double
    end
    methods (Abstract)
        getManifoldSize(this);
    end
    methods
        function p = pdf(this, xa)
            % Use nearest neighbor interpolation by default
            [~, indices] = getClosestPoint(this, xa);
            p = this.gridValues(indices)';
        end
        
        function grid = getGrid(this)
            % Overload if .grid should stay empty
            grid = this.grid;
        end

        function gridPoints = getGridPoint(this, indices)
            % To avoid passing all points if only one or few are needed.
            % Overload if .grid should stay empty
            arguments
                this (1,1) AbstractGridDistribution
                indices double {mustBeInteger} = []
            end
            if isempty(indices)
                gridPoints = this.grid;
            else
                gridPoints = this.grid(:,indices);
            end
        end
        
        function [points, indices] = getClosestPoint(this, xa)
            % Overload if class does not have .grid
            allDistances = angularError(reshape(this.grid,this.dim,1,[]),xa);
            if this.dim>1 % Combine into single dimension for multidimensional case
                % This is good for both hypertori and hyperspheres (the
                % ordering is the same as with the orthodromic distance)
                % Not for hyperhemispheres (is overloaded there)
                allDistances = vecnorm(allDistances,2,1);
            end
            [~,indices] = min(allDistances,[],3);
            points = this.getGridPoint(indices);
        end
        
        function gd = multiply(this, other)
            arguments
                this AbstractGridDistribution
                other AbstractGridDistribution
            end
            assert(this.enforcePdfNonnegative == other.enforcePdfNonnegative);
            gd = this;
            gd.gridValues = gd.gridValues.*other.gridValues;
            gd = gd.normalize(warnUnnorm=false);
        end
        
        function gd = normalize(this, opt)
            arguments
                this AbstractGridDistribution
                opt.tol (1,1) double = 1e-4
                opt.warnUnnorm (1,1) logical = true
            end
            int = this.integral();
            gd = this; % Asssign here so we can use return in the else branch
            if any(this.gridValues < 0)
                warning('Normalization:negative', 'There are negative values. This usualy points to a user error.');
            elseif abs(int) < 1e-200 % Tolerance has to be that low to avoid unnecessary errors on multiply
                error('Normalization:almostZero', 'Sum of grid vallues is too close to zero, this usually points to a user error.');
            elseif abs(int-1) > opt.tol
                if opt.warnUnnorm
                    warning('Normalization:notNormalized', 'Grid values apparently do not belong to normalized density. Normalizing...');
                end
            else % Density is normalized. Do not change anything. We need this because we also want to to normalize in all other cases
                return
            end
            gd.gridValues = gd.gridValues/int;
        end
        
        function int = integral(this, l, r)
            arguments
                this (1,1) AbstractGridDistribution
                l double = []
                r double = []
            end
            assert(isempty(l)==isempty(r), 'Either both limits or only none should be empty.');
            if isempty(l)&&~strcmp(this.gridType,'sh_grid') % Currently only grid that surely requires different integral
                int = this.getManifoldSize * mean(this.gridValues);
            elseif isempty(l)
                warning('Using numerical integral for this type of grid.');
                int = this.integralNumerical();
            else
                warning('Grid:CanOnlyIntegrateWithLimitsNumerically','Integration with limits currently only supported numerically.')
                int = this.integralNumerical(l, r);
            end
        end
                
    end
end

