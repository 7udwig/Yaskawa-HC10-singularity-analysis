% dh_transform - Denavit-Hartenberg Transformation Matrix
%
% Computes the 4x4 homogeneous transformation matrix A_i for a single
% robot joint using the standard Denavit-Hartenberg (DH) convention.
%
% The DH transformation is defined as a sequence of two translations
% and two rotations:
%   A = Tz(d) * Rz(theta) * Tx(a) * Rx(alpha)
%
% This yields the matrix:
%   | cos(theta)  -sin(theta)*cos(alpha)   sin(theta)*sin(alpha)  a*cos(theta) |
%   | sin(theta)   cos(theta)*cos(alpha)  -cos(theta)*sin(alpha)  a*sin(theta) |
%   |     0             sin(alpha)              cos(alpha)              d       |
%   |     0                 0                       0                  1       |
%
% INPUTS:
%   d      - Link offset along the previous z-axis [mm]
%   thetar - Joint angle about the previous z-axis [rad]
%   a      - Link length along the x-axis (perpendicular distance
%            between z-axes) [mm]
%   alpha  - Link twist angle about the x-axis [rad]
%
% OUTPUT:
%   DHT    - 4x4 homogeneous transformation matrix from frame i-1 to frame i
%
% USAGE EXAMPLE:
%   T = dh_transform(275, deg2rad(0), 0, -pi/2);
%
% See also: forward_kinematics, singularity_analysis

function [DHT] = dh_transform(d, thetar, a, alpha)

    DHT = [cos(thetar), -sin(thetar)*cos(alpha),  sin(thetar)*sin(alpha), a*cos(thetar);
           sin(thetar),  cos(thetar)*cos(alpha), -cos(thetar)*sin(alpha), a*sin(thetar);
           0,            sin(alpha),              cos(alpha),             d            ;
           0,            0,                       0,                      1            ];

end
