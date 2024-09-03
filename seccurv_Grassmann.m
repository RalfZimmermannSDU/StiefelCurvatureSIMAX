function [seccurv] = seccurv_Grassmann(B1,B2)
% Compute the sectional curvature associated with the plane
% spanned by X= [0;B1], Y=[0;B2] in T_[I;0]Gr(n,p)
%   
%INPUT
% B1,B2 : "coordinates" of tangent vectors on 
%         Grassmann "mathfrak m", linear independent.
%OUTPUT
% seccurv : sectional curvature of span(X,Y)
%

% 1 step normalize tangents
B1 = B1/norm(B1,'fro');
B2 = B2-trace(B2'*B1)*B1;
%checkorth = trace(B1'*B2)
B2 = B2/norm(B2, 'fro');

B12 = (B1'*B2);

M1221 = B12*B12';
M1122 = (B1'*B1)*(B2'*B2);
M1212 = B12*B12;

seccurv = trace(M1221) + trace(M1122) - 2*trace(M1212);

% for comparison purposes: alternative way to comute the curvature
%compare = 0.5*(norm(-B1'*B2 + B2'*B1,'fro')^2 + norm(-B1*B2' + B2*B1','fro')^2)
end