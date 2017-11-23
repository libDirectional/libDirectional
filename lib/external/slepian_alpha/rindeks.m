function x=rindeks(M,i)
% x=RINDEKS(M,i)
%
% Returns the ROW i of a matrix M.
% Returns the ith 1-st dimension of a 3D matrix.
%
% See also INDEKS, KINDEKS.
%
% Last modified by fjsimons-at-alum.mit.edu, 12/03/2009

x=M(i(:),:,:);
