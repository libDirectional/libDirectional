classdef CustomToroidalDistribution < AbstractToroidalDistribution & CustomHypertoroidalDistribution
    % Toroidal distribution with custom pdf.
    
    methods
        function this = CustomToroidalDistribution(f_)
            % It is the user's responsibility to ensure that f is a valid
            % circular density, i.e., 2pi-periodic, nonnegative and
            % normalized.
            this@CustomHypertoroidalDistribution(f_,2);
        end
    end
    
end

