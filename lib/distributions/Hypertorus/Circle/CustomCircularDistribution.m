classdef CustomCircularDistribution < AbstractCircularDistribution & CustomHypertoroidalDistribution
    % Circular distribution with custom pdf.
    
    methods
        function this = CustomCircularDistribution(f_)
            % It is the user's responsibility to ensure that f is a valid
            % circular density, i.e., 2pi-periodic, nonnegative and
            % normalized.
            this@CustomHypertoroidalDistribution(f_,1);
        end
    end
    
end

