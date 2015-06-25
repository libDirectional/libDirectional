% Calculates product of complex numbers q and w
% Parameters:
%   q (2 x 1 column vector)
%       first complex number
%   w (2 x 1 column vector)
%       second complex number
% Returns:
%   r (2 x 1 column vector)
%       product q*w

function r = complexMultiplication(q, w)
    assert(all(size(q) == [2,1]));
    assert(all(size(w) == [2,1]));
    
    r = [q(1,:) .* w(1,:) - q(2,:) .* w(2,:); q(1,:) .* w(2,:) + q(2,:) .* w(1,:)];
end
