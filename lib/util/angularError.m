%% Calculates the angular error between alpha and beta
function e = angularError(alpha, beta)
    arguments
        alpha double {mustBeNonempty}
        beta double {mustBeNonempty}
    end
    alpha = mod(alpha,2*pi);
    beta = mod(beta,2*pi);
    diff = abs(alpha-beta);
    e = min(diff,2*pi-diff);
end