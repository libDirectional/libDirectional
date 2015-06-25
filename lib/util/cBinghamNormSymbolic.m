function [c, X] = cBinghamNormSymbolic(D, debug)
    % This function is a symbolic calculation of the complex Bingham
    % normalization constant.
    %
    % Returns symbolic normalization constant and vector of eigenvalues.
    
    if nargin == 0
        clear all;
        clc;
        D = 3;
        debug = true;
    elseif nargin == 1
        debug = false;
    end
    
    %% Define variables
    % Vector of eigenvalues
    X = sym('x', [D, 1]);
    
    % Vector of prefactos
    B = sym(ones(D, 1));
    
    %% Define a vector of differences
    dims = 1:D;
    
    for d = dims
        for dd = dims(dims~=d)
            B(d) = B(d) * (X(d) - X(dd));
        end
    end
    
    if debug
        pretty(B);
    end
    
    %% Define prefactors
    B = 1 ./ B;
    
    if debug
        pretty(B);
    end
    
    %% Normalization constant
    p = sym(pi*ones(1, 1));
    c = 2 * p^D * sum(B .* exp(X));
    
    if debug
        pretty(c);
    end
end

