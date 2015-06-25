classdef AbstractToroidalFilter < handle
    % Abstract base class for filters on the torus
    
    properties
    end
    
	methods (Abstract)
        setState(this, state)
        est = getEstimate(this)
    end
    
end

