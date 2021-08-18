classdef HemisphericalGridDistribution < HyperhemisphericalGridDistribution
    methods
        function hgd = toFullSphere(this)
            % Convert hemisphere to full sphere (grid is symmetric but
            % values may not ensure symmetry due to numerical issues
            grid_ = [this.grid,-this.grid];
            gridValues_ = 0.5*[this.gridValues;this.gridValues];
            hgd = SphericalGridDistribution(grid_,gridValues_);
        end
        
        function h = plotInterpolated(this, useHarmonics)
            arguments
                this HyperhemisphericalGridDistribution
                useHarmonics (1,1) logical = true
            end
            hdgd = this.toFullSphere();
            hhgdInterp = CustomHemisphericalDistribution(@(x)2*hdgd.pdf(x,useHarmonics));
            warnStruct = warning('off', 'PDF:UseInterpolated');
            h = hhgdInterp.plot;
            warning(warnStruct);
        end
    end
    methods (Static)
        function sgd = fromDistribution(dist, noOfGridPoints, gridType)
            arguments
                dist AbstractDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set_symmetric'
            end
            hgd = HyperhemisphericalGridDistribution.fromDistribution(dist, noOfGridPoints, gridType);
            sgd = HemisphericalGridDistribution(hgd.grid, hgd.gridValues, hgd.enforcePdfNonnegative);
            sgd.gridType = hgd.gridType;
        end
        function sgd = fromFunction(fun, noOfGridPoints, gridType)
            arguments
                fun function_handle
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set_symm'
            end
            % Always use adjusted version of eq_point_set since only this is suited to the hemihypersphere
            hgd = HyperhemisphericalGridDistribution.fromFunction(fun, noOfGridPoints, 3, gridType);
            sgd = HemisphericalGridDistribution(hgd.grid, hgd.gridValues, hgd.enforcePdfNonnegative);
            sgd.gridType = hgd.gridType;
        end
    end
end