classdef CustomCircularDistribution < AbstractCircularDistribution
    % Circular distribution with custom pdf.
    
    properties
        f
    end
    
    methods
        function this = CustomCircularDistribution(f_)
            % It is the user's responsibility to ensure that f is a valid
            % circular density, i.e., 2pi-periodic, nonnegative and
            % normalized.
            assert(isa(f_, 'function_handle'));
            this.f = f_;
        end
        
        function p = pdf(this, xa)
            p = arrayfun(this.f, xa);
        end
    end
    
end

