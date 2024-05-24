clear; close all;

%% define rotation matrices
theta = sym('theta', [3, 1]);

Rx = [1 0 0; ...
      0 cos(theta(1)) -sin(theta(1)); ...
      0 sin(theta(1)) cos(theta(1))];

Ry = [cos(theta(2)) 0 sin(theta(2)); ...
      0 1 0; ...
      -sin(theta(2)) 0 cos(theta(2))];

Rz = [cos(theta(3)) -sin(theta(3)) 0; ...
      sin(theta(3)) cos(theta(3)) 0; ...
      0 0 1];

%% rotation combination
R = Rx*Ry*Rz;

matlabFunction(R, 'Vars', {theta}, 'File', 'eul2rotm_func');

