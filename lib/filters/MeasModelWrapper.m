classdef MeasModelWrapper < AdditiveNoiseMeasurementModel
    % Helper class for using the nonlinear filtering toolbox, converts a
    % measurement function f into an additive noise measurement model.
    
    properties
        f function_handle
    end
    
    methods
        function this = MeasModelWrapper(f_)
            this.f = f_;
        end
        
        function measurements = measurementEquation(this, states, ~)
            for i=1:size(states,2)
                measurements(:,i) = this.f(states(:,i)); %#ok<AGROW>
            end
        end
    end
    
end

