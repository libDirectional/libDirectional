classdef ToroidalUniformDistribution < AbstractToroidalDistribution & HypertoroidalUniformDistribution
    % Uniform distribution on the torus
    
    methods
        function this = ToroidalUniformDistribution()
            % Constructor
            this@HypertoroidalUniformDistribution(2);
        end
        
        function tu = shift(~,~)
            tu = ToroidalUniformDistribution;
        end
    end
end

