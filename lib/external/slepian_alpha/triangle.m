function t=triangle(l1,l2,l3)
% t=TRIANGLE(l1,l2,l3)
%
% Checks the triangle inequality, Dahlen and Tromp (1998), Eq. (C.186)
% If t is FALSE then the Wigner3j symbols are zero.
%
% INPUT:
%
% l1, l2, l3    Same-sized vectors of degrees
%
% Last modified by fjsimons-at-alum.mit.edu, 11/19/2010

l1=l1(:);
l2=l2(:);
l3=l3(:);

if ~all([length(l1) length(l2) length(l3)]==length(l1)); 
  error('All inputs must have same length') ; 
end

% Don't do the double operand
t=abs(l2-l3)<=l1 & l1<=(l2+l3);
