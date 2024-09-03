
%
%
%
% This script conducts the experiments in Section 4.3 
% and produces Figures 3 and 4 of the SIMAX paper
%
% HIGH CURVATURE MEANS LOW-RANK: ON THE SECTIONAL
% CURVATURE OF GRASSMANN AND STIEFEL MANIFOLDS AND
% THE UNDERLYING MATRIX TRACE INEQUALITIES
%
% by RALF ZIMMERMANN AND JAKOB STOYE
%
%
%



clear;
close all;
% ------------
%
% The three matrices from the DDVV paper
% "DDVV-type inequality for skew-symmetric matrices and
% Simons-type inequality for Riemannian submersions"
% by Jianquan Ge
% ------------

A1 = @(s) [[0,s,0,0];[-s,0,0,0];[0,0,0,s];[0,0,-s,0]];
A2 = @(s) [[0,0,s,0];[0,0,0,-s];[-s,0,0,0];[0,s,0,0]];
A3 = @(s) [[0,0,0,s];[0,0,s,0];[0,-s,0,0];[-s,0,0,0]];

% for rank tests
B1 = @(u,v) [[0,1,0,0];[-1,0,0,0];[0,0,0,u];[0,0,-v,0]];
B2 = @(u,v) [[1,0,0,0];[0,1,0,0];[0,0,u,0];[0,0,0,v]];





% test dimensions p = 4, n = 8
%
res = 51;   % plot resolution
unitI = linspace(0,1,res);

curv_AvsB_cplot_A1A2  = zeros(res,1);   % canonical 
curv_AvsB_eplot_A1A2  = zeros(res,1);   % Euclidean 
curv_AvsB_cplot_A3A2  = zeros(res,1);   % canonical 
curv_AvsB_eplot_A3A2  = zeros(res,1);   % Euclidean
curv_rankB_cplot = zeros(res,res);      % canonical
curv_rankB_eplot = zeros(res,res);      % Euclidean

for k = 1:res
    u = unitI(k);
    % data for Figure 4
    curv_AvsB_cplot_A1A2(k) = seccurv_Stiefel_canon( u*A1(1), (1-u)*B1(0,0), u*A2(1), (1-u)*B2(0,0));
    curv_AvsB_eplot_A1A2(k) = seccurv_Stiefel_euclid(u*A1(1), (1-u)*B1(0,0), u*A2(1), (1-u)*B2(0,0));
    curv_AvsB_cplot_A3A2(k) = seccurv_Stiefel_canon( u*A3(1), (1-u)*B1(0,0), u*A2(1), (1-u)*B2(0,0));
    curv_AvsB_eplot_A3A2(k) = seccurv_Stiefel_euclid(u*A3(1), (1-u)*B1(0,0), u*A2(1), (1-u)*B2(0,0));

    for j=1:res
        v = unitI(j);
        %data for Figure 3
        curv_rankB_cplot(k,j) = seccurv_Stiefel_canon(zeros(4), B2(u,v), zeros(4),  B1(u,v));
        curv_rankB_eplot(k,j) = seccurv_Stiefel_euclid(zeros(4), B2(u,v), zeros(4),  B1(u,v));
    end
end

figure;  % produce Figure 3 of the paper
subplot(1,2,1)
[unitX, unitY] = meshgrid(unitI, unitI);
surf(unitX,unitY, curv_rankB_cplot,'edgecolor','none')
axis([0 1 0 1 0 1.25])
ax = gca;
ax.FontSize = 20; 

title('Sectional curvature (canonical) for B-block depending on u and v')
xlabel('$u$', 'interpreter','latex','FontSize',14)
ylabel('$v$', 'interpreter','latex','FontSize',14)
zlabel('$K(0,B_1(u,v), 0, B_2(u,v))$', 'interpreter','latex','FontSize',24)
hold on
subplot(1,2,2)
[unitX, unitY] = meshgrid(unitI, unitI);
surf(unitX,unitY, curv_rankB_eplot,'edgecolor','none')
axis([0 1 0 1 0 1.25])
ax = gca;
ax.FontSize = 20; 

title('Sectional curvature (Euclid) for B-block depending on u and v')
xlabel('$u$', 'interpreter','latex','FontSize',14)
ylabel('$v$', 'interpreter','latex','FontSize',14)
zlabel('$K(0,B_1(u,v), 0, B_2(u,v))$', 'interpreter','latex','FontSize',24)
%-----------------------------------

figure; % produce Figure 4 of the paper
subplot(1,2,1)
plot(unitI,curv_AvsB_cplot_A1A2, 'k',...
     unitI,curv_AvsB_eplot_A1A2, 'k:',LineWidth=2)
ylim([0 1.5])
yticks(linspace(0,1.5,7))
ax = gca;
ax.FontSize = 20; 
title('Sectional curvature for increasing A-block')
xlabel('$u$', 'interpreter','latex','FontSize',24)
ylabel('$K(uA_1,(1-u)B_1, uA_2, (1-u)B_2)$', 'interpreter','latex','FontSize',24)
legend("Stiefel (can.)", "Stiefel (Euclid)")

subplot(1,2,2)
plot(unitI,curv_AvsB_cplot_A3A2, 'k',...
     unitI,curv_AvsB_eplot_A3A2, 'k:',LineWidth=2)
ylim([0 1.5])
yticks(linspace(0,1.5,7))
ax = gca;
ax.FontSize = 20; 
title('Sectional curvature for increasing A-block')
xlabel('$u$', 'interpreter','latex','FontSize',24)
ylabel('$K(uA_3,(1-u)B_1, uA_2, (1-u)B_2)$', 'interpreter','latex','FontSize',24)
legend("Stiefel (can.)", "Stiefel (Euclid)")