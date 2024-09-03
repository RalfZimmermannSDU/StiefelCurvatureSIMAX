function [seccurv] = seccurv_SOn(X,Y)
% Compute the sectional curvature associated with the plane
% spanned by X= [A1, -B1'; B1, C1], Y=[A2,-B2';B2, C2] in T_[I]SOn
%   
%INPUT
% X,Y : "coordinates" of tangent vectors on 
%         Grassmann "mathfrak m", linear independent.
%OUTPUT
% seccurv : sectional curvature of span(X,Y)
%

% 1 step normalize tangents
X = X/norm(X,'fro');
Y = Y-trace(X'*Y)*X;
%checkorth = trace(X'*Y)
Y = Y/norm(Y, 'fro');

Lie_XY = X*Y-Y*X;
seccurv = 0.5*trace(Lie_XY'*Lie_XY);
end