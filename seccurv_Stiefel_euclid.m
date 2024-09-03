function [seccurv, X,Y] = seccurv_Stiefel_euclid(A1, B1, A2, B2)
% Compute the sectional curvature associated with the plane
%                |A1|      |A2|
% spanned by X = |B1|, Y = |B2| in T_[I]St(n,p)
%       EUCLIDEAN METRIC
%INPUT
% A1,A2,B1,B2 : "coordinates" of tangent vectors on 
%                Stiefel "mathfrak m", linear independent.
%                A1,A2 skew
%OUTPUT
% seccurv : sectional curvature of span(X,Y)
% X,Y     : orthogonalized X,Y

% dimensions
p  = size(A1,1);
np = size(B1,1);  % n minus p

% 1 step
%   orthonormalize tangents w.r.t. Euclidean metric

normX = sqrt(trace(A1'*A1) + trace(B1'*B1));
A1 = A1/normX;
B1 = B1/normX;
% keep orthogonal component of Y
d = trace(A1'*A2) + trace(B1'*B2);
A2 = A2 - d*A1;
B2 = B2 - d*B1;
normY = sqrt(trace(A2'*A2) + trace(B2'*B2));
A2 = A2/normY;
B2 = B2/normY;

LieA1A2  = A1*A2  - A2*A1;
LieB1TB2 = B1'*B2 - B2'*B1;
LieB2B1T = B2*B1' - B1*B2';

seccurv = (1/4)*norm(LieA1A2 + LieB1TB2,'fro')^2 ...
        +       norm(B1*A2-B2*A1, 'fro')^2 ...
        +       trace(B1*(B2'*B2)*B1') - trace((B1'*B2)*(B2'*B1));

% ONB of section
X = [A1;B1];
Y = [A2;B2];
end