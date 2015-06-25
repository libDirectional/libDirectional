classdef AbstractCircularFilter < handle
    % Abstract base class for filters on the circle (S1, SO(2))
    
    properties
    end
    
    methods (Abstract)
        setState(this, state)
        est = getEstimate(this)
    end
    
end

