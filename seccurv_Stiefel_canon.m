function [seccurv, X,Y] = seccurv_Stiefel_canon(A1, B1, A2, B2)
% Compute the sectional curvature associated with the plane
%                |A1  -B1'|      |A2   -B2'|
% spanned by X = |B1    0 |, Y = |B2     0 | in T_[I]St(n,p)
%
% CANONICAL METRIC
%
%
%INPUT
% A1,A2,B1,B2 : "coordinates" of tangent vectors on 
%                Stiefel "mathfrak m", linear independent.
%                A1,A2 skew
%OUTPUT
% seccurv : sectional curvature of span(X,Y)
%

% dimensions
p  = size(A1,1);
np = size(B1,1);  % n minus p

% 1 step
%   orthonormalize tangents w.r.t. canonical metric

normX = sqrt(0.5*trace(A1'*A1) + trace(B1'*B1));
A1 = A1/normX;
B1 = B1/normX;
% keep orthogonal component of Y
d = 0.5*trace(A1'*A2) + trace(B1'*B2);
A2 = A2 - d*A1;
B2 = B2 - d*B1;
normY = sqrt(0.5*trace(A2'*A2) + trace(B2'*B2));
A2 = A2/normY;
B2 = B2/normY;

LieA1A2  = A1*A2  - A2*A1;
LieB1TB2 = B1'*B2 - B2'*B1;
LieB2B1T = B2*B1' - B1*B2';

seccurv = (1/8)*norm(LieA1A2 - LieB1TB2,'fro')^2 ...
        + (1/4)*norm(B1*A2-B2*A1, 'fro')^2 ...
        + (1/2)*norm(LieB2B1T, 'fro')^2;

% construct matrices that span ONB of tangent plane section
X = [[A1,-B1'];[B1,zeros(size(B1,1))]];
Y = [[A2,-B2'];[B2,zeros(size(B2,1))]];
end