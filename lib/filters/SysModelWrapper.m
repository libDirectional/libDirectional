classdef SysModelWrapper < AdditiveNoiseSystemModel
    % Helper class for using the nonlinear filtering toolbox, converts a
    % system function f into an additive noise system model.    
    
    properties
       f function_handle
    end 
    
    methods
        function this = SysModelWrapper(f_)
            this.f=f_;
        end
        
        function predictedStates = systemEquation(this, states)
            for i=1:size(states,2)
                predictedStates(:,i) = this.f(states(:,i)); %#ok<AGROW>
            end
        end
    end
end