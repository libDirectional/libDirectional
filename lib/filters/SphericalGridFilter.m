classdef SphericalGridFilter < AbstractHypersphericalFilter & HypersphericalGridFilter
    
    methods
        function this = SphericalGridFilter(noOfCoefficients, gridType)
            % Constructor
            %
            % Parameters:
            %   noOfCoefficients (integer > 0)
            %       number of grid values to use
            %	gridType (char)
            %       grid type to use for the SphericalGridDistribution
            arguments
                noOfCoefficients {mustBeInteger,mustBePositive}
                gridType char = 'eq_point_set'
            end
            % Only allow eq_point_set because large errors are made in the
            % prediction step for the spherical harmonics-based grid.
            this@HypersphericalGridFilter(noOfCoefficients, 3, gridType);
            this.gd = SphericalGridDistribution.fromDistribution(HypersphericalUniformDistribution(3),...
                noOfCoefficients, gridType);
        end
        
        function setState(this, gd_)
            % Sets the current system state
            %
            % Parameters:
            %   gd_ (AbstractHypersphericalDistribution)
            %       new state
            % This function is only overloaded to verify the class of the
            % argument
            setState@HypersphericalGridFilter(this, gd_);
            assert(isa(this.gd,'SphericalGridDistribution')); % Verify afterward to allow for conversion in setState of the superclass function
        end
    end
    
end
