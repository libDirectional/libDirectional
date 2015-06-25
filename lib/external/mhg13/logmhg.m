% function s=logmhg(M,alpha,p,q,x,y)
%
% computes the logarithm of the 
% truncated hypergeometric function pFq ^alpha(p;q;x;y) 
% of one OR two matrix arguments
% with the sum over |kappa|<=M 
% p,q,x,y are arrays, so hg(6,9,[3 4],[5 6 7],[1,2],[8 9]) is 
% 2F3^9([3 4];[5 6 7];[1 2];[8,9]) with the sum over all partitions
% |kappa|<=6
% y may be omitted
%
% Plamen Koev
% Department of Mathematics
% Massachusetts Institute of Technology
%
% Written November 2004
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for details

error('No mex file found for this platform');                    
