classdef S2CondS2GridDistribution < SdCondSdGridDistribution
    % In this class, a conditional distribution on the sphere is
    % described by grid values. For the conditional distribution f(a|b), we
    % allow both a and by to vary. Thus, it is a function of the Cartesian
    % product of two spheres. It obviously should only integrate to 1 for a
    % fixed b. To provide a grid for the Cartesian product of two spheres, we
    % generate the Cartesian product of the grids of the individual spheres.
    methods
        function this = S2CondS2GridDistribution(grid_, gridValues_)
            % Conditional distribution! First dim conditioned on second
            % dim.
            % Provide grid on the sphere, Cartesian product
            % will be grid on S2 x S2
            assert(size(grid_,1)==3);
            this = this@SdCondSdGridDistribution(grid_, gridValues_);
        end     
    end
    methods (Static)
        function s2s2 = fromFunction(fun, noGridPoints, funDoesCartesianProduct, gridType)
            arguments
                fun function_handle
                noGridPoints (1,1) {mustBeInteger}
                % State if function does Cartesian product itself. Use this
                % if keeping one argument constant should be done by the
                % function, e.g., because it can be realized more
                % efficiently.
                funDoesCartesianProduct logical = false
                gridType char = 'eq_point_set'
            end
            sdsd = SdCondSdGridDistribution.fromFunction(fun, noGridPoints, funDoesCartesianProduct, gridType, 6);
            s2s2 = S2CondS2GridDistribution(sdsd);
        end
    end
end

