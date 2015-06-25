% function [s,c]=mhgi(MAX,alpha,p,q,n,x)
%
% computes the truncated hypergeometric function pFq ^alpha(p;q;x*I_n)
% The sum is only over |kappa|<=MAX(1) and kappa(i)<=MAX(2)
% if size(MAX)=[1,1], then MAX(2) is by default equal to MAX(1) 
% p and q are arrays, so hgi(30,9,[3 4],[5 6 7],3,0.5) is 
% 2F3^9([3 4],[5 6 7];0.5*I_3) summed over all kappa with |kappa|<30
% Uses the formula 
% J_lambda(x I_n)=x^|lambda| \prod_{(i,j)\in \lambda (n-(i-1)+alpha(j-1))
%
% When the value of mhgi is desired for multiple values of x one 
% can input x as a vector. Then mhgi returns a vector with the values of 
% mhgi for every component of x.
% 
% s is a polynomial in x. The output vector c contains the coefficients of
% that polynomial so that s=c(1)+x*c(2)+...+x^M*c(M+1), meant componentwise
% when x is a vector. Equivalently s=polyval(c(end:-1:1),x)
%
% Plamen Koev
% Department of Mathematics
% Massachusetts Institute of Technology
%
% May 2004
% Copyright (c) 2004 Plamen Koev. See COPYRIGHT.TXT for details
%
% Reference: Plamen Koev and Alan Edelman, The Efficient Evaluation of the
% Hypergeometric Function of a Matrix Argument, Math. Comp. 75 (2006),
% 833-846.

error('No mex file found for this platform');
