function grad_log_c = cBinghamGradLogNormSymbolic(c, X, debug)
    % Symbolic calculation of first derivative of logarithm of normalization
    % constant of the complex Bingham distribution.
    
    if nargin == 2
        debug = false;
    end
    
    %% Gradient of logarithm of normalization constant (for ML estimation)
    grad_log_c = jacobian(log(c), X);
    
    if debug
        pretty(grad_log_c);
    end
end

