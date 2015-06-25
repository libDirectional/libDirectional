function samples = discretesample(p, n)
    % Obtains samples from discrete distribution given by p
    %
    % Parameters:
    %   p (row or column vector)
    %       propability for each discrete possibility (need not be
    %       normalized)
    %   n (integer > 0)
    %       number of samples
    p=p(:);
    assert(all(p>=0));
    assert(n>0);
    
    P = cumsum(p);
    total = P(end);
    assert(total>0);
    [~,samples] = histc(rand(1,n), [0; P]/total);
end