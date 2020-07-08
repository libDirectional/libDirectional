classdef SphericalGridDistribution < AbstractHypersphericalDistribution & HypersphericalGridDistribution
    
    methods
        function this = SphericalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_)
            % Constructor
            arguments % Use to set default value
                grid_ double {mustBeLessThanOrEqual(grid_,1),mustBeGreaterThanOrEqual(grid_,-1)}
                gridValues_ (:,1) double {mustBeNonnegative}
                enforcePdfNonnegative_ logical = true
            end
            this = this@HypersphericalGridDistribution(grid_, gridValues_, enforcePdfNonnegative_);
            % Validate it has the right dimension
            assert(this.dim==3);
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
                    degree = (-6+sqrt(36-8*(4-noOfGridPoints)))/4; % Größeres Ergebnis der abc-Formel
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
            sgd = SphericalGridDistribution(grid, gridValues, enforcePdfNonnegative);
        end
    end
end

