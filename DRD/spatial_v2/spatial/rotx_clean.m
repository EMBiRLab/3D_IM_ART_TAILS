function  X = rotx_clean( theta, threshold )
% if theta is very close to pi, pi/2, 0,
% directly return the integer matrix
% to reduce rounding error in reading urdfs

% rotx  spatial coordinate transform (X-axis rotation).
% rotx(theta)  calculates the coordinate transform matrix from A to B
% coordinates for spatial motion vectors, where coordinate frame B is
% rotated by an angle theta (radians) relative to frame A about their
% common X axis.

if nargin < 2
    threshold = 1e-8;
end

if abs(angdiff(theta, pi)) <= threshold
    c = -1;
    s = 0;
elseif abs(angdiff(theta, pi/2)) <= threshold
    c = 0;
    s = 1;
elseif abs(angdiff(theta, 0)) <= threshold
    c = 1;
    s = 0;
elseif abs(angdiff(theta, -pi/2)) <= threshold
    c = 0;
    s = -1;
else
    c = cos(theta);
    s = sin(theta);
end

X = [ 1  0  0  0  0  0 ;
      0  c  s  0  0  0 ;
      0 -s  c  0  0  0 ;
      0  0  0  1  0  0 ;
      0  0  0  0  c  s ;
      0  0  0  0 -s  c
    ];
