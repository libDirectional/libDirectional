function points_x = eq_point_set_symm(dim, N, deleteHalf, symmetryType, varargin)
%EQ_POINT_SET Center points of regions of EQ partition, in Cartesian coordinates
%
%Syntax
% points_x = eq_point_set(dim,N,options);
%
%Description
% POINTS_X = EQ_POINT_SET(dim,N) does the following:
% 1) uses the recursive zonal equal area sphere partitioning algorithm to
% partition S^dim (the unit sphere in dim+1 dimensional space) into N regions
% of equal area and small diameter, and
% 2) sets POINTS_X to be an array of size (dim+1 by N), containing the center
% points of each region.
% Each column of POINTS_X represents a point of S^dim, in Cartesian coordinates.
%
% The arguments dim and N must be positive integers.
%
% POINTS_X = EQ_POINT_SET(dim,N,'offset','extra') uses experimental extra offsets
% for S^2 and S^3 to try to minimize energy.
%
% POINTS_X = EQ_POINT_SET(dim,N,extra_offset) uses experimental extra offsets if
% extra_offset is true or non-zero.
%
%Notes
% Each region is defined as a product of intervals in spherical polar
% coordinates. The center point of a region is defined via the center points
% of each interval, with the exception of spherical caps and their descendants,
% where the center point is defined using the center of the spherical cap.
%
% If dim > 3, extra offsets are not used.
% For more details on options, see help partition_options.
%
%Examples
% > points_x = eq_point_set(2,4)
% points_x =
%          0    0.0000   -0.0000    0.0000
%          0    1.0000   -1.0000         0
%     1.0000    0.0000    0.0000   -1.0000
%
% > size(points_x)
% ans =
%      3     4
%
%See also
% PARTITION_OPTIONS, EQ_POINT_SET_POLAR, EQ_REGIONS, S2X

% Adapted version by Florian Pfaff (pfaff@kit.edu) for libDirectional
% 2020-03-28
% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Function changed name from s2x to polar2cart
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

%
% Check number of arguments
%
arguments
    dim (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(dim,2)}
    N (1,1) double {mustBeInteger(N),mustBeGreaterThanOrEqual(N,2)}
    deleteHalf (1,1) logical = false
    symmetryType char = 'mirror' % Can be 'mirror' or 'plane'
end
arguments (Repeating)
    varargin
end
if deleteHalf
    points_x = polar2cart(eq_point_set_polar_plane_symm(dim, N, true, varargin{:}));
elseif strcmp(symmetryType,'mirror')
    points_x = polar2cart(eq_point_set_polar_plane_symm(dim, N, true, varargin{:}));
    points_x = [points_x, -points_x];
elseif strcmp(symmetryType,'plane')
    points_x = polar2cart(eq_point_set_polar_plane_symm(dim, N, false, varargin{:}));
elseif strcmp(symmetryType,'asymm')
    points_x = polar2cart(eq_point_set_polar(dim, N, varargin{:}));
else
    error('Unknown sysmmetriy type');
end
