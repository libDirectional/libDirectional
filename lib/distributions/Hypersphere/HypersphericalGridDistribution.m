classdef HypersphericalGridDistribution < AbstractHypersphericalDistribution & AbstractGridDistribution
	properties
        grid double {mustBeLessThanOrEqual(grid,1),mustBeGreaterThanOrEqual(grid,-1)}
	end
    methods
        function this = HypersphericalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            assert(size(grid_,2)==size(gridValues_,1));
            if size(grid_,1)>size(grid_,2)
                warning('Dimension is higher than number of grid points. Verify that is really intended.');
            end
            this.dim = size(grid_,1);
            this.grid = grid_;
            this.gridValues = gridValues_;
            this.enforcePdfNonnegative = enforcePdfNonnegative_;
            % Check if normalized. If not: Normalize!
            this = this.normalize;
        end      
        
        function mu = meanDirection(this)
            mu = sum(this.grid.*this.gridValues',2);
            if norm(mu)<1e-8
                warning('Density may not have actually have a mean direction because formula yields a point very close to the origin.')
            end
            mu = mu/norm(mu);
        end
        
        function f = normalize(this)
            tol = 1e-2;
            f = normalize@AbstractGridDistribution(this,tol);
        end
                
        function hgd = multiply(this, other)
            arguments
                this HypersphericalGridDistribution
                other HypersphericalGridDistribution
            end
            assert(isequal(this.grid,other.grid), 'Multiply:IncompatibleGrid','Can only multiply for equal grids.');
            hgd = multiply@AbstractGridDistribution(this, other);
        end   
        
        function p = pdf(this, xa, useHarmonics)
            arguments
                this HypersphericalGridDistribution
                xa double
                useHarmonics (1,1) logical = true
            end
            warning('PDF:UseInterpolated', 'pdf is not defined. Using interpolation with spherical harmonics.')
            if useHarmonics
                if this.enforcePdfNonnegative
                    transformation = 'sqrt';
                else
                    transformation = 'identity';
                end
                shd = SphericalHarmonicsDistributionComplex.fromGrid(this.gridValues, this.grid, transformation);
                p = shd.pdf(xa);
            else
                warning('Interpolating the pdf is not very efficient, but it is good enough for visualization purposes.');
                [~,maxIndex]=max(this.grid'*xa);
                p = this.gridValues(maxIndex)';
            end
        end
        
        function h = plot(this)
            if this.dim == 3
                AbstractHypersphericalDistribution.plotSphere;
            end
            hold on
            hdd = HypersphericalDiracDistribution(this.grid, this.gridValues');
            h = hdd.plot;
            hold off
        end
        
        function h = plotInterpolated(this, useHarmonics)
            arguments
                this HypersphericalGridDistribution
                useHarmonics (1,1) logical = true
            end
            if this.dim~=3
                error('Can currently only plot for S2 sphere.')
            end
            if useHarmonics
                if this.enforcePdfNonnegative
                    transformation = 'sqrt';
                else
                    transformation = 'identity';
                end
                shd = SphericalHarmonicsDistributionComplex.fromGrid(this.gridValues, this.grid, transformation);
                h = shd.plot;
            else
                chd = CustomHypersphericalDistribution(@(x)this.pdf(x,false),3);
                chd.plot;
            end
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
        function sgd = fromFunction(fun, noOfGridPoints, dim, gridType)
            arguments
                fun function_handle
                noOfGridPoints {mustBeInteger,mustBePositive}
                dim {mustBeInteger,mustBeGreaterThanOrEqual(dim,2)}
                gridType char = 'eq_point_set'
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
            sgd = HypersphericalGridDistribution(grid, gridValues);
        end
    end
end

