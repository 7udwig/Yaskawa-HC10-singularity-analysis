% forward_kinematics - Forward Kinematics for Yaskawa HC10 DTP Robot
%
% Computes the TCP (Tool Center Point) pose as a 4x4 homogeneous
% transformation matrix by chaining the six individual DH transformation
% matrices for the Yaskawa HC10 DTP robot.
%
% The DH parameters used here were determined empirically to match
% the MotoSim simulation environment (absolute coordinate system):
%
%   Joint | d [mm] | theta offset | a [mm] | alpha [rad]
%   ------|--------|--------------|--------|------------
%     1   |  275   |      0       |    0   |   -pi/2
%     2   |    0   |    -90°      |  700   |    pi
%     3   |    0   |      0       |    0   |   -pi/2
%     4   | -500   |      0       |    0   |    pi/2
%     5   | -162   |      0       |    0   |   -pi/2
%     6   | -170   |      0       |    0   |    pi
%
% NOTE: Joint 2 requires a -90° offset to match the robot's zero position.
%       This was verified against Yaskawa MotoSim.
%
% INPUTS:
%   A1..A6 - Joint angles [rad] for joints 1 through 6
%            (apply any offsets BEFORE passing to this function)
%
% OUTPUT:
%   TCP_Pose - 4x4 homogeneous transformation matrix of the TCP
%              The translational part (column 4, rows 1-3) gives the
%              TCP position in [mm]. The rotational part (rows 1-3,
%              cols 1-3) gives the orientation.
%
% USAGE EXAMPLE:
%   angles_deg = [0, -90, 0, 0, 0, 0];   % joint angles with offsets
%   angles_rad = deg2rad(angles_deg);
%   T = forward_kinematics(angles_rad(1), angles_rad(2), angles_rad(3), ...
%                          angles_rad(4), angles_rad(5), angles_rad(6));
%   px = T(1,4);  % TCP x-position [mm]
%   py = T(2,4);  % TCP y-position [mm]
%   pz = T(3,4);  % TCP z-position [mm]
%
% See also: dh_transform, singularity_analysis

function TCP_Pose = forward_kinematics(A1, A2, A3, A4, A5, A6)

    % --- DH Parameters for Yaskawa HC10 DTP (verified against MotoSim) ---
    %
    % d:     link offsets along z-axis [mm]
    % a:     link lengths along x-axis [mm]
    % alpha: link twist angles [rad]
    %
    % Note: theta offsets (-90° on joint 2) must be applied BEFORE
    %       calling this function (see singularity_analysis.m for the convention used).

    d     = [  275;    0;    0;  -500;  -162;  -170];
    a     = [    0;  700;    0;     0;     0;     0 ];
    alpha = [-pi/2;   pi; -pi/2; pi/2; -pi/2;  pi  ];

    % Compute individual DH transformation matrices for each joint
    T1 = dh_transform(d(1), A1, a(1), alpha(1));
    T2 = dh_transform(d(2), A2, a(2), alpha(2));
    T3 = dh_transform(d(3), A3, a(3), alpha(3));
    T4 = dh_transform(d(4), A4, a(4), alpha(4));
    T5 = dh_transform(d(5), A5, a(5), alpha(5));
    T6 = dh_transform(d(6), A6, a(6), alpha(6));

    % Chain all transformations: base -> joint1 -> ... -> TCP
    TCP_Pose = T1 * T2 * T3 * T4 * T5 * T6;

end
