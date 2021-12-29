classdef HypercylindricalStateSpaceSubdivisionDistribution < StateSpaceSubdivisionDistribution & AbstractHypercylindricalDistribution
    methods
        function this = HypercylindricalStateSpaceSubdivisionDistribution(gd_, linDistributions)
            arguments
                gd_ (1,1) AbstractGridDistribution
                linDistributions (:,1) AbstractLinearDistribution
            end
            this@StateSpaceSubdivisionDistribution(gd_, linDistributions);
        end
        
        function h = plot(this, interpolate)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionDistribution
                interpolate (1,1) logical = false
            end
            if interpolate % Cannot call plotInterpolated directly because only superclass methods with same name can be called
                h = plot@AbstractHypercylindricalDistribution(this);
            else
                h = plot@StateSpaceSubdivisionDistribution(this);
            end
        end
        
        function h = plotInterpolated(this)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionDistribution
            end
            h = plot(this, true);
        end
        
        function m = mode(this)
            arguments
                this (1,1) HypercylindricalStateSpaceSubdivisionDistribution
            end
            m = mode@StateSpaceSubdivisionDistribution(this);
        end
    end
    
    methods (Static)
        function hcrbd = fromDistribution(dist, noOfGridPoints, gridType)
            arguments
                dist AbstractHypercylindricalDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'CartesianProd'
            end
            hcrbd = HypercylindricalStateSpaceSubdivisionDistribution.fromFunction(@(x)dist.pdf(x),...
                noOfGridPoints, dist.linD, dist.boundD, gridType);
        end
        function hcrbd = fromFunction(fun, noOfGridPoints, dimLin, dimBound, gridType, intRange)
            arguments
                fun (1,1) function_handle
                noOfGridPoints (1,1) {mustBeInteger,mustBePositive}
                dimLin (1,1) {mustBeInteger,mustBePositive}
                dimBound (1,1) {mustBeInteger,mustBePositive} = 1 %#ok<INUSA>
                gridType char = 'CartesianProd' %#ok<INUSA>
                intRange (1,:) double = [-inf,inf]
            end
            assert(nargin(fun) == 1, 'Need to be given in the format used for .pdf in densities.');
            assert(dimLin == 1, 'Currently, bounded dimension must be 1.');
            
            gd = FIGDistribution.fromDistribution(CircularUniformDistribution, noOfGridPoints);
            grid = gd.getGrid();
            % Preallocate
            cds = repmat(CustomLinearDistribution(@(x)x,1),[1,noOfGridPoints]);
            
            for i=1:noOfGridPoints
                funCurr = @(y)reshape(fun([grid(i).*ones(1,numel(y)); y(:)']),size(y));
                % Obtain grid value via integral
                gd.gridValues(i) = integral(funCurr,intRange(1),intRange(2));
                % Original function divided by grid value is linear
                cds(i) = CustomLinearDistribution(@(x)funCurr(x)/gd.gridValues(i),1);
            end
            
            hcrbd = HypercylindricalStateSpaceSubdivisionDistribution(gd,cds);
        end
    end
end