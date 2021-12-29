classdef HypersphericalGridDistribution < AbstractHypersphereSubsetGridDistribution & AbstractHypersphericalDistribution
	methods
        function this = HypersphericalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_, gridType)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
                gridType char = 'unknown'
            end
            this@AbstractHypersphereSubsetGridDistribution(grid_, gridValues_, enforcePdfNonnegative_);
            this.gridType = gridType;
        end      
        
        function mu = meanDirection(this)
            mu = sum(this.grid.*this.gridValues',2);
            if norm(mu)<1e-8
                warning('Density may not actually have a mean direction because formula yields a point very close to the origin.')
            end
            mu = mu/norm(mu);
        end

        function p = pdf(this, xa)
            arguments
                this HypersphericalGridDistribution
                xa double
            end
            warning('PDF:UseInterpolated','Interpolating the pdf with constant values in each region is not very efficient, but it is good enough for visualization purposes.');
            [~,maxIndex]=max(this.grid'*xa);
            p = this.gridValues(maxIndex)';
        end
        
        function h = plot(this)
            if this.dim == 3
                AbstractHypersphericalDistribution.plotSphere;
            end
            hold on
            hdd = HypersphericalDiracDistribution(this.grid, this.gridValues'/sum(this.gridValues));
            h = hdd.plot;
            hold off
        end
        
        function h = plotInterpolated(this)
            arguments
                this HypersphericalGridDistribution
            end
            chd = CustomHypersphericalDistribution(@(x)this.pdf(x,false),this.dim);
            h = chd.plot;
        end
        
        function hhgd = symmetrize(this)
            assert(isequal(this.grid(:,2),-this.grid(:,2+size(this.grid,2)/2)),'Symmetrize:AsymmetricGrid',...
                'Can only use symmetrize for asymetric grids. Use eq_point_set_symm when calling fromDistribution or fromFunction.');
            gridValuesHalf = 0.5 * (this.gridValues(1:size(this.grid,2)/2) + this.gridValues(size(this.grid,2)/2+1:end));
            hhgd = HypersphericalGridDistribution(this.grid, [gridValuesHalf; gridValuesHalf]);
        end
        
        function hhgd = toHemisphere(this, tol)
            arguments
                this HypersphericalGridDistribution
                tol (1,1) double = 1e-10
            end
            assert(isequal(this.grid(:,2),-this.grid(:,2+size(this.grid,2)/2)),'ToHemisphere:AsymmetricGrid',...
                'Can only use symmetrize for asymetric grids. Use eq_point_set_symm when calling fromDistribution or fromFunction.');
            if abs(this.gridValues(2)-this.gridValues(2+size(this.grid,2)/2))<tol
                gridValuesHemisphere = 2*(this.gridValues(1:size(this.grid,2)/2));
            else
                warning('ToHemisphere:AsymmetricDensity','Density appears to be asymmetric. Not converting to hemispherical one.');
                % 2 cancels out with 0.5
                gridValuesHemisphere = (this.gridValues(1:size(this.grid,2)/2) + this.gridValues(size(this.grid,2)/2+1:end));
            end
            hhgd = HypersphericalGridDistribution(this.grid, 2*gridValuesHemisphere);
        end
    end
    
    methods (Static)
        function sgd = fromDistribution(dist, noOfGridPoints, gridType)
            arguments
                dist AbstractHypersphericalDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set'
            end
            sgd = HypersphericalGridDistribution.fromFunction(@(x)dist.pdf(x), noOfGridPoints, dist.dim, gridType);
        end
        function sgd = fromFunction(fun, noOfGridPoints, dim, gridType, enforcePdfNonnegative)
            arguments
                fun function_handle
                noOfGridPoints {mustBeInteger,mustBePositive}
                dim {mustBeInteger,mustBeGreaterThanOrEqual(dim,2)}
                gridType char = 'eq_point_set'
                enforcePdfNonnegative logical = true
            end
            switch gridType
                case 'eq_point_set'
                    grid = eq_point_set(dim-1, noOfGridPoints);
                case {'eq_point_set_symm','eq_point_set_symmetric'}
                    grid = eq_point_set_symm(dim-1, noOfGridPoints, false, 'mirror');
                case 'eq_point_set_symm_plane'
                    grid = eq_point_set_symm(dim-1, 2*noOfGridPoints, false, 'plane');
                otherwise
                    error('Grid scheme not recognized');
            end
            gridValues = fun(grid)';
            sgd = HypersphericalGridDistribution(grid, gridValues, enforcePdfNonnegative, gridType);
        end
    end
end

