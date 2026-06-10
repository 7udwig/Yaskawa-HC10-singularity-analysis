% singularity_analysis - Geometric Jacobian and Singularity Analysis for Yaskawa HC10 DTP
%
% This script reads joint angle data from a CSV/Excel file, computes the
% forward kinematics via DH transformations, assembles the geometric
% Jacobian matrix for each configuration, and evaluates its determinant
% as a measure of kinematic manipulability.
%
% A determinant value close to zero indicates a kinematic singularity,
% meaning the robot loses one or more degrees of freedom in that pose.
%
% WORKFLOW:
%   1. Load joint angle data from file (one row per timestep, 6 columns)
%   2. For each timestep:
%       a. Apply joint angle offsets (DH convention)
%       b. Compute forward kinematics (DH chain)
%       c. Extract z-axes and positions from transformation matrices
%       d. Build the geometric Jacobian column by column
%       e. Compute and store the determinant
%   3. Export results (TCP positions + determinant) as CSV files
%
% INPUT FILE FORMAT:
%   - CSV or Excel file with N rows and 6 columns
%   - Each column corresponds to one joint angle [degrees]: q1 ... q6
%   - No header row expected (raw data only)
%   - Default filename: 'joint_angles.csv'
%     (can be overridden by setting INPUT_FILE below)
%
% OUTPUT FILES (written to current directory):
%   - determinant.csv   : Jacobian determinant per timestep
%   - position_x.csv    : TCP x-position per timestep [mm]
%   - position_y.csv    : TCP y-position per timestep [mm]
%   - position_z.csv    : TCP z-position per timestep [mm]
%
% DEPENDENCIES:
%   dh_transform.m, forward_kinematics.m
%
% USAGE:
%   1. Place your joint angle CSV file in the same folder as this script.
%   2. Set INPUT_FILE to your filename (or leave as default).
%   3. Run the script. Results are saved as CSV files.
%
% EXAMPLE JOINT ANGLE OFFSETS (DH convention, applied before calculation):
%   theta_DH = [theta1, theta2 - 90, theta3, theta4, theta5, theta6]
%
% See also: dh_transform, forward_kinematics

clc;
clear;

% =========================================================================
% --- CONFIGURATION -------------------------------------------------------
% =========================================================================

INPUT_FILE = 'joint_angles.csv';   % Input file (CSV or .xlsx)
                                   % 6 columns: q1..q6 in degrees, no header

% =========================================================================
% --- LOAD DATA -----------------------------------------------------------
% =========================================================================

fprintf('Loading joint angle data from: %s\n', INPUT_FILE);

% Support both CSV and Excel input
[~, ~, ext] = fileparts(INPUT_FILE);
if strcmpi(ext, '.csv')
    M = readmatrix(INPUT_FILE);
elseif strcmpi(ext, '.xlsx') || strcmpi(ext, '.xls')
    M = xlsread(INPUT_FILE);   % xlsread kept for legacy Excel compatibility
else
    error('Unsupported file format: %s. Use .csv or .xlsx.', ext);
end

% Validate that the file has at least 6 columns
if size(M, 2) < 6
    error('Input file must have at least 6 columns (one per joint angle).');
end

numRows = size(M, 1);
fprintf('Loaded %d data rows.\n', numRows);

% =========================================================================
% --- DH PARAMETERS (Yaskawa HC10 DTP, verified against MotoSim) ---------
% =========================================================================
%
% The first coordinate system is placed at the intersection of axes 1 and 2
% (Yaskawa convention: d1 = 275 mm). The -90° offset on joint 2 aligns the
% DH model with the robot's absolute zero position.
%
%   Joint | d [mm] | theta offset | a [mm] | alpha [rad]
%   ------|--------|--------------|--------|------------
%     1   |  275   |      0       |    0   |   -pi/2
%     2   |    0   |    -90°      |  700   |    pi
%     3   |    0   |      0       |    0   |   -pi/2
%     4   | -500   |      0       |    0   |    pi/2
%     5   | -162   |      0       |    0   |   -pi/2
%     6   | -170   |      0       |    0   |    pi

d     = [  275;    0;    0;  -500;  -162;  -170];
a     = [    0;  700;    0;     0;     0;     0 ];
alpha = [-pi/2;   pi; -pi/2; pi/2; -pi/2;  pi  ];

% =========================================================================
% --- PRE-ALLOCATE OUTPUT ARRAYS ------------------------------------------
% =========================================================================

Det    = zeros(numRows, 1);   % Jacobian determinant
pos_x  = zeros(numRows, 1);   % TCP x-position [mm]
pos_y  = zeros(numRows, 1);   % TCP y-position [mm]
pos_z  = zeros(numRows, 1);   % TCP z-position [mm]

% =========================================================================
% --- MAIN LOOP -----------------------------------------------------------
% =========================================================================

fprintf('Processing configurations...\n');

for i = 1:numRows

    % --- Read joint angles for this timestep [degrees] -------------------
    theta1 = M(i, 1);
    theta2 = M(i, 2);
    theta3 = M(i, 3);
    theta4 = M(i, 4);
    theta5 = M(i, 5);
    theta6 = M(i, 6);

    % --- Apply DH angle offsets and convert to radians -------------------
    % Joint 2 requires -90° offset to match the MotoSim zero position.
    thetaDeg = [theta1; theta2 - 90; theta3; theta4; theta5; theta6];
    A = deg2rad(thetaDeg);

    % --- Compute individual DH transformation matrices -------------------
    T1 = dh_transform(d(1), A(1), a(1), alpha(1));
    T2 = dh_transform(d(2), A(2), a(2), alpha(2));
    T3 = dh_transform(d(3), A(3), a(3), alpha(3));
    T4 = dh_transform(d(4), A(4), a(4), alpha(4));
    T5 = dh_transform(d(5), A(5), a(5), alpha(5));
    T6 = dh_transform(d(6), A(6), a(6), alpha(6));

    % --- Cumulative transformation matrices (base -> joint i) ------------
    TT1 = T1;
    TT2 = T1 * T2;
    TT3 = T1 * T2 * T3;
    TT4 = T1 * T2 * T3 * T4;
    TT5 = T1 * T2 * T3 * T4 * T5;
    TT6 = T1 * T2 * T3 * T4 * T5 * T6;

    % --- Extract TCP position from final transformation matrix -----------
    px = TT6(1, 4);
    py = TT6(2, 4);
    pz = TT6(3, 4);

    % --- Extract origin positions of each frame --------------------------
    % Column 4 of each cumulative matrix gives the frame's origin in base coords.
    t0 = [0; 0; 0];               % Base frame origin
    t1 = TT1(1:3, 4);
    t2 = TT2(1:3, 4);
    t3 = TT3(1:3, 4);
    t4 = TT4(1:3, 4);
    t5 = TT5(1:3, 4);
    t6 = TT6(1:3, 4);             % TCP position

    % --- Extract z-axis unit vectors for each frame ----------------------
    % Column 3 of each cumulative matrix is the z-axis direction in base coords.
    z0 = [0; 0; 1];               % Base frame z-axis
    z1 = TT1(1:3, 3);
    z2 = TT2(1:3, 3);
    z3 = TT3(1:3, 3);
    z4 = TT4(1:3, 3);
    z5 = TT5(1:3, 3);

    % --- Build the geometric Jacobian column by column -------------------
    % For a revolute joint i, the Jacobian column is:
    %   J_i = [ z_{i-1} x (p_TCP - p_{i-1}) ]
    %         [          z_{i-1}             ]
    % where x denotes the cross product.

    J1 = [cross(z0, t6 - t0); z0];
    J2 = [cross(z1, t6 - t1); z1];
    J3 = [cross(z2, t6 - t2); z2];
    J4 = [cross(z3, t6 - t3); z3];
    J5 = [cross(z4, t6 - t4); z4];
    J6 = [cross(z5, t6 - t5); z5];

    Jacobian = [J1, J2, J3, J4, J5, J6];   % 6x6 geometric Jacobian

    % --- Compute determinant as singularity measure ----------------------
    % det(J) -> 0 indicates a kinematic singularity (loss of mobility).
    determinant = det(Jacobian);

    % --- Store results ---------------------------------------------------
    Det(i)   = determinant;
    pos_x(i) = px;
    pos_y(i) = py;
    pos_z(i) = pz;

    % Progress feedback every 1000 rows
    if mod(i, 1000) == 0
        fprintf('  Processed %d / %d rows...\n', i, numRows);
    end

end

fprintf('Done. Writing output files...\n');

% =========================================================================
% --- EXPORT RESULTS ------------------------------------------------------
% =========================================================================

writematrix(Det,   'determinant.csv');
writematrix(pos_x, 'position_x.csv');
writematrix(pos_y, 'position_y.csv');
writematrix(pos_z, 'position_z.csv');

fprintf('Results saved:\n');
fprintf('  determinant.csv   - Jacobian determinant per timestep\n');
fprintf('  position_x.csv    - TCP x-position [mm]\n');
fprintf('  position_y.csv    - TCP y-position [mm]\n');
fprintf('  position_z.csv    - TCP z-position [mm]\n');

% =========================================================================
% --- OPTIONAL: QUICK PLOT ------------------------------------------------
% =========================================================================

figure;
subplot(2,1,1);
plot(pos_x, 'r'); hold on;
plot(pos_y, 'g');
plot(pos_z, 'b');
xlabel('Sample index');
ylabel('Position [mm]');
title('TCP Position over Time');
legend('X', 'Y', 'Z');
grid on;

subplot(2,1,2);
plot(Det, 'k');
xlabel('Sample index');
ylabel('det(J)');
title('Jacobian Determinant (singularity indicator)');
yline(0, 'r--', 'Singularity');
grid on;
