classdef S2CondS2GridDistribution < SdCondSdGridDistribution
    
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

