classdef HyperhemisphericalGridDistribution < AbstractHypersphereSubsetGridDistribution & AbstractHyperhemisphericalDistribution
    methods
        function this = HyperhemisphericalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            assert(all(grid_(end,:)>=0), 'Always using upper hemisphere (along last dimension).');
            this@AbstractHypersphereSubsetGridDistribution(grid_, gridValues_, enforcePdfNonnegative_);
        end      
        
        function mu = meanDirection(this)
            warning('For hyperhemispheres, this function yields the mode and not the mean.');
            % If we took the mean, it would be biased toward [0;...;0;1] 
            % because the lower half is considered inexistant.
            [~,indexMax] = max(this.gridValues);
            mu = this.grid(:,indexMax);
        end
                
        function hgd = toFullSphere(this)
            % Convert hemisphere to full sphere (grid is symmetric but
            % values may not ensure symmetry due to numerical issues
            grid_ = [this.grid,-this.grid];
            gridValues_ = 0.5*[this.gridValues;this.gridValues];
            hgd = HypersphericalGridDistribution(grid_,gridValues_);
        end

        function h = plot(this)
            hdd = HypersphericalDiracDistribution(this.grid, this.gridValues');
            h = hdd.plot;
        end
        
        function h = plotInterpolated(this)
            arguments
                this HyperhemisphericalGridDistribution
            end
            hdgd = this.toFullSphere;
            hhgdInterp = CustomHyperhemisphericalDistribution(@(x)2*hdgd.pdf(x),3);
            warnStruct = warning('off', 'PDF:UseInterpolated');
            h = hhgdInterp.plot;
            warning(warnStruct);
        end
        
        function h = plotFullSphereInterpolated(this)
            if this.dim~=3
                error('Can currently only plot for hemisphere of S2 sphere.')
            end
            hgd = this.toFullSphere;
            shd = SphericalHarmonicsDistributionComplex.fromGrid(hgd.gridValues, hgd.grid, 'identity');
            chhd = CustomHypersphericalDistribution.fromDistribution(shd);
            h = chhd.plot;
        end
        
        function [points, indices] = getClosestPoint(this, xa)
            % Consider symmetry by using minimum. We can always assume it's
            % multidimensional so we always use vecnorm. For more
            % explanation see the definition in AbstractGridDistribution.
            allDistances = min(...
                vecnorm(angularError(reshape(this.grid,this.dim,1,[]),xa),2,1),...
                vecnorm(angularError(reshape(-this.grid,this.dim,1,[]),xa),2,1));
            [~,indices] = min(allDistances,[],3);
            points = this.getGridPoint(indices);
        end
    end
    
    methods (Static)
        function sgd = fromDistribution(dist, noOfGridPoints, gridType)
            arguments
                dist AbstractDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set_symmetric'
            end
            if isa(dist,'AbstractHyperhemisphericalDistribution')
                fun = @(x)dist.pdf(x);
            elseif isa(dist,'WatsonDistribution') || isa(dist,'VMFDistribution') && dist.mu(end)==0 || isa(dist,'BinghamDistribution') ||...
                    isa(dist,'HypersphericalMixture')... % Allow for mixtures of two von Mises--Fisher
                    && numel(dist.dists)==2 && all(dist.w==0.5) && isequal(dist.dists{1}.mu,-dist.dists{2}.mu)
                fun = @(x)2*dist.pdf(x);
            elseif isa(dist,'HypersphericalGridDistribution')
                error('Converting a HypersperhicalGridDistribution to a HypersphericalGridDistribution is not supported');
            elseif isa(dist,'AbstractHypersphericalDistribution')
                warning('fromDistribution:asymmetricOnHypersphere',...
                    'Approximating a hyperspherical distribution on a hemisphere. The density may not be symmetric. Double check if this is intentional.');
                fun = @(x)2*dist.pdf(x);
            else
                error('Distribution currently not supported.');
            end
            sgd = HyperhemisphericalGridDistribution.fromFunction(fun, noOfGridPoints, dist.dim, gridType);
        end
        function sgd = fromFunction(fun, noOfGridPoints, dim, gridType)
            arguments
                fun function_handle
                noOfGridPoints {mustBeInteger,mustBePositive}
                dim {mustBeInteger,mustBeGreaterThanOrEqual(dim,2)}
                gridType char = 'eq_point_set_symm'
            end
            % Always use adjusted version of eq_point_set since only this is suited to the hemihypersphere
            switch gridType
                case {'eq_point_set', 'eq_point_set_symm', 'eq_point_set_symmetric'}
                    grid = eq_point_set_symm(dim-1, 2*noOfGridPoints, true);
                otherwise
                    error('Grid scheme not recognized');
            end
            gridValues = fun(grid)';
            sgd = HyperhemisphericalGridDistribution(grid, gridValues);
        end
    end
end

