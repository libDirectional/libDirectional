classdef AbstractHypersphericalFilter < handle
    % Abstract base class for filters on the hypersphere (without antipodal
    % symmetry). Use AbstractAxialFilter for antipodally symmetric
    % problems!
    
    properties
        dim %dimension (dim=2 -> circle/complex numbers, dim=4 -> quaternions, etc.)
    end
        
    methods (Abstract)
        setState(this, state)
        est = getEstimate(this)
        mean=getEstimateMean(this)
    end
       
end
