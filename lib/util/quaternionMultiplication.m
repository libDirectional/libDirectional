% Calculates product of quaternions q and w
% Parameters:
%   q (4 x 1 column vector)
%       first quaternion
%   w (4 x 1 column vector)
%       second quaternion
% Returns:
%   r (4 x 1 column vector)
%       product q*w
function r = quaternionMultiplication(q, w)    
    assert(all(size(q) == [4 1]));
    assert(all(size(w) == [4 1]));
    
    % Formula based on cross product
    % Hart, J. C.; Francis, G. K. & Kauffman, L. H. 
    % Visualizing Quaternion Rotation ACM Transactions on Graphics, 
    % ACM, 1994, 13, 256-276
    
    q0 = q(1,:);
    w0 = w(1,:);
    qvec = q(2:4,:);
    wvec = w(2:4,:);
    r = [q0 .* w0 - dot(qvec, wvec); q0 .* wvec + w0 .* qvec + cross(qvec, wvec)];
end