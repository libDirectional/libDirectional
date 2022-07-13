classdef AbstractHypersphereSubsetGridDistribution < AbstractGridDistribution & AbstractHypersphereSubsetDistribution
    methods
        function this = AbstractHypersphereSubsetGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            assert(size(grid_,2)==size(gridValues_,1));
            this.dim = size(grid_,1);
            this.grid = grid_;
            this.gridValues = gridValues_;
            this.enforcePdfNonnegative = enforcePdfNonnegative_;
            % Check if normalized. If not: Normalize!
            this = this.normalize;
            this.gridType = 'unknown';
        end      
        
        function mu = meanDirection(this)
            warning('For hyperhemispheres, this function yields the mode and not the mean.');
            % If we took the mean, it would be biased toward [0;...;0;1] 
            % because the lower half is considered inexistant.
            [~,indexMax] = max(this.gridValues);
            mu = this.grid(:,indexMax);
        end

        function C = moment(this)
            C = this.grid.*(this.gridValues'/sum(this.gridValues))*this.grid';
        end
        
        function f = normalize(this, opt)
            arguments
                this (1,1) AbstractHypersphereSubsetGridDistribution
                opt.tol (1,1) double = 1e-2
                opt.warnUnnorm (1,1) logical = true
            end
            f = normalize@AbstractGridDistribution(this,tol=opt.tol,warnUnnorm=opt.warnUnnorm);
        end
        
        function int = integral(this)
            arguments
                this (1,1) AbstractHypersphereSubsetGridDistribution
            end
            int = integral@AbstractGridDistribution(this);
        end
        
        function hgd = multiply(this, other)
            arguments
                this AbstractHypersphereSubsetGridDistribution
                other AbstractHypersphereSubsetGridDistribution
            end
            assert(isequal(this.grid,other.grid), 'Multiply:IncompatibleGrid','Can only multiply for equal grids.');
            hgd = multiply@AbstractGridDistribution(this, other);
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

