function grad_c_divided_by_c = cBinghamGradNormDividedByNormSymbolic(c, X, debug)
    % Symbolic calculation of the first derivative of the normalization constant
    % divided by the normalization constant itself.
    
    if nargin == 2
        debug = false;
    end
    
    %% Gradient
    grad_c_divided_by_c = jacobian(c, X);
    grad_c_divided_by_c = grad_c_divided_by_c ./ c;
    
    if debug
        pretty(grad_c_divided_by_c);
    end
end

