%
%
%
% This script conducts the experiments in Section 4.2 
% and produces Figure 2 of the SIMAX paper
%
% HIGH CURVATURE MEANS LOW-RANK: ON THE SECTIONAL
% CURVATURE OF GRASSMANN AND STIEFEL MANIFOLDS AND
% THE UNDERLYING MATRIX TRACE INEQUALITIES
%
% by RALF ZIMMERMANN AND JAKOB STOYE
% 
%
% the output data is written to file
%


runs    = input("How many runs for random experiments? ")
% setting in the paper: runs=100;
max_dim = input("Maximum dimension p in experiments? ")  
% setting in the paper: max_dim=1000;
start_dim = 2;
% initialize curvature arrays
curv_rank_Stiefel_cplot  = zeros(max_dim,1);
curv_rank_Stiefel_eplot  = zeros(max_dim,1);
curv_rank_Grassmann_plot = zeros(max_dim,1);
curv_rank_SOn_plot       = zeros(max_dim,1);

for p = start_dim:max_dim
    for k = 1:runs
        % create random tangents
        X = rand(2*p);
        X = 0.5*(X-X');
        Y = rand(2*p);
        Y = 0.5*(Y'-Y);
        Ax = X(1:p,1:p);
        Bx = X(p+1:2*p,1:p);
        Ay = Y(1:p,1:p);
        By = Y(p+1:2*p,1:p);

        % compute curvatures
        curv_rank_Stiefel_cplot(p)  = curv_rank_Stiefel_cplot(p) + ...
                                       seccurv_Stiefel_canon(Ax,Bx,Ay,By);
        curv_rank_Stiefel_eplot(p)  = curv_rank_Stiefel_eplot(p) + ...
                                       seccurv_Stiefel_euclid(Ax,Bx,Ay,By); 
        curv_rank_Grassmann_plot(p)= curv_rank_Grassmann_plot(p) +...
                                       seccurv_Grassmann(Bx,By);
        curv_rank_SOn_plot(p)      = curv_rank_SOn_plot(p) + ...
                                       seccurv_SOn(X,Y);
    end
    % take average over runs
    curv_rank_Stiefel_cplot(p)  = curv_rank_Stiefel_cplot(p)/runs;
    curv_rank_Stiefel_eplot(p)  = curv_rank_Stiefel_eplot(p)/runs;
    curv_rank_Grassmann_plot(p) = curv_rank_Grassmann_plot(p)/runs;
    curv_rank_SOn_plot(p)       = curv_rank_SOn_plot(p)/runs;
end


figure;

semilogy(start_dim:max_dim, curv_rank_Grassmann_plot(start_dim:max_dim), 'b-.', ...
     start_dim:max_dim, curv_rank_Stiefel_cplot(start_dim:max_dim), 'k',...
     start_dim:max_dim, curv_rank_Stiefel_eplot(start_dim:max_dim), 'k:',...
     start_dim:max_dim, curv_rank_SOn_plot(start_dim:max_dim), 'r--', LineWidth=2)
ax = gca;
ax.FontSize = 24; 
legend("Grassmann", "Stiefel (can.)", "Stiefel (Euclid)", "SO(n)")
title('Sectional curvature versus increasing dimension')
xlabel('$p$', 'interpreter','latex','FontSize',24)
ylabel('Sectional curvature $K$', 'interpreter','latex','FontSize',24)

% write data to disk
BigRankTest.Grassmann = curv_rank_Grassmann_plot;
BigRankTest.Stiefelc  = curv_rank_Stiefel_cplot;
BigRankTest.Stiefele  = curv_rank_Stiefel_eplot;
BigRankTest.SOn  = curv_rank_SOn_plot;
eval(['save BigRankTest_maxdim', num2str(max_dim), '.mat BigRankTest']);

