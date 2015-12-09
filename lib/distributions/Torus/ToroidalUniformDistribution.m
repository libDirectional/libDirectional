classdef ToroidalUniformDistribution < AbstractToroidalDistribution & HypertoroidalUniformDistribution
    % Uniform distribution on the torus
    
    methods
        function this = ToroidalUniformDistribution()
            % Constructor
            this@HypertoroidalUniformDistribution(2);
        end
        
        function cu = shift(~,~)
            cu = ToroidalUniformDistribution;
        end
    end
end

