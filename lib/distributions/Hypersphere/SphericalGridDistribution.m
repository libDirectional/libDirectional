classdef SphericalGridDistribution < AbstractHypersphericalDistribution & HypersphericalGridDistribution
    
    methods
        function this = SphericalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_, gridType)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
                gridType char = 'unknown'
            end
            this = this@HypersphericalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_, gridType);
            % Validate it has the right dimension
            assert(this.dim==3);
        end
        
        function f = normalize(this,opt)
            arguments
                this (1,1) SphericalGridDistribution
                opt.tol (1,1) double = 1e-2
                opt.warnUnnorm (1,1) logical = true
            end
            if ~strcmp(this.gridType,'sh_grid')
                f = normalize@AbstractGridDistribution(this,tol=opt.tol,warnUnnorm=opt.warnUnnorm);
            else
                warning('SphericalGridDistribution:CannotNormalizeShGrid','Cannot properly normalize for sh_grid');
                f = normalize@AbstractGridDistribution(this,tol=opt.tol,warnUnnorm=opt.warnUnnorm);
            end
        end
        
        function h = plotInterpolated(this, useHarmonics)
            arguments
                this HypersphericalGridDistribution
                useHarmonics (1,1) logical = true
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
        
        function p = pdf(this, xa, useHarmonics)
            arguments
                this HypersphericalGridDistribution
                xa double
                useHarmonics (1,1) logical = true
            end
            if useHarmonics
                warning('PDF:UseInterpolated', 'pdf is not defined. Using interpolation with spherical harmonics.')
                if this.enforcePdfNonnegative
                    transformation = 'sqrt';
                else
                    transformation = 'identity';
                end
                if ~strcmp(this.gridType,'sh_grid')
                    shd = SphericalHarmonicsDistributionComplex.fromGrid(this.gridValues, this.grid, transformation);
                else
                    warnStruct = warning('off','Normalization:notNormalized');
                    warning('PDF:NeedNormalizationShGrid','Need to normalize for pdf because it is unnormalized since sh_grid is used.');
                    shd = SphericalHarmonicsDistributionComplex.fromGrid(this.gridValues, this.grid, transformation);
                    warning(warnStruct);
                end
                
                p = shd.pdf(xa);
            else
                warning('PDF:UseInterpolated','Interpolating the pdf with constant values in each region is not very efficient, but it is good enough for visualization purposes.');
                [~,maxIndex]=max(this.grid'*xa);
                p = this.gridValues(maxIndex)';
            end
        end
    end
    
    methods (Static)
        function sgd = fromDistribution(dist, noOfGridPoints, gridType, enforcePdfNonnegative)
            arguments
                dist AbstractHypersphericalDistribution
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set'
                enforcePdfNonnegative logical = true
            end
            assert(dist.dim==3,'fromDistribution:WrongDimension',...
                'Distribution must have dimension 3 to approximate using SphericalGridDistribution.');
            sgd = SphericalGridDistribution.fromFunction(@(x)dist.pdf(x), noOfGridPoints, gridType, enforcePdfNonnegative);
        end
        function sgd = fromFunction(fun, noOfGridPoints, gridType, enforcePdfNonnegative)
            arguments
                fun function_handle
                noOfGridPoints {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set'
                enforcePdfNonnegative logical = true
            end
            switch gridType
                case 'eq_point_set'
                    grid = eq_point_set(2, noOfGridPoints);
                case 'sh_grid'
                    warning('Transformation:notEq_Point_set',...
                        'Not using eq_point_set. This may lead to problems in the normalization (and filters based thereon should not be used because the transition may not be valid).');
                    degree = (-6+sqrt(36-8*(4-noOfGridPoints)))/4; % Larger solution of the quadratic equation
                    assert(degree==round(degree), 'Number of coefficients not supported for this type of grid.');
                    lat = linspace(0, 2*pi, 2*degree+2);
                    lon = linspace(pi/2, -pi/2, degree+2);
                    [latMesh, lonMesh] = meshgrid(lat, lon);
                    [x, y, z] = sph2cart(latMesh(:)', lonMesh(:)', 1);
                    grid = [x;y;z];
                otherwise
                    error('Grid scheme not recognized');
            end
            gridValues = fun(grid)';
            sgd = SphericalGridDistribution(grid, gridValues, enforcePdfNonnegative, gridType);
        end
    end
end

