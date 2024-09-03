%
%
%
% This script conducts the experiments in Section 4.1 
% and produces Figure 1 of the SIMAX paper
%
% HIGH CURVATURE MEANS LOW-RANK: ON THE SECTIONAL
% CURVATURE OF GRASSMANN AND STIEFEL MANIFOLDS AND
% THE UNDERLYING MATRIX TRACE INEQUALITIES
%
% by RALF ZIMMERMANN AND JAKOB STOYE
% 
%
%
%
%
n = 10;          % set dimension
A1 = zeros(n);   % initialize matrix blocks
A2 = zeros(n);
B1 = zeros(n);
B2 = zeros(n);


B1(1,2) = 1;     % set start matrices as displayed 
B2(1,1) = 1;     % in Section 4.1 of the paper

res = 1001;      % resolution for plots
unitI = linspace(0,1,res);

% initialize arrays to record the curvature
curv_rank_Stiefel_cplot  = zeros(9*res,1); % canonical Stiefel
curv_rank_Stiefel_eplot  = zeros(9*res,1); % Euclidean Stiefel
curv_rank_Stiefel_eTplot = zeros(9*res,1); % Euclidean Stiefel (transposed B)
curv_rank_Grassmann_plot = zeros(9*res,1); % Grassmann curvature
curv_rank_SOn_plot       = zeros(9*res,1); % SO(n) curvature
for k = 1:9*res
    %
    if k<= res
        % fill the u2-entries
        B1(2,1) = -unitI(k);
        B2(2,2) =  unitI(k);
        % compute curvature
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=2*res
        % keep u2, fill u3
        B2(3,3) =  unitI((k-res));
        B1(3,4) =  unitI((k-res));
        % compute curvature
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=3*res
        % keep u2,u3, fill u4
        B1(4,3) = -unitI((k-2*res));
        B2(4,4) =  unitI((k-2*res));
        % compute curvatures
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=4*res
        % keep u2,...,u4, fill u5
        B1(5,6) =  unitI((k-3*res));
        B2(5,5) =  unitI((k-3*res));
        % compute curvatures
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=5*res
        % keep u2,...,u5, fill u6
        B1(6,5) = -unitI((k-4*res));
        B2(6,6) =  unitI((k-4*res));
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=6*res
        % keep u2,...,u6, fill u7
        B1(7,8) =  unitI((k-5*res));
        B2(7,7) =  unitI((k-5*res));
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=7*res
        % keep u2,...,u7, fill u8
        B1(8,7) = -unitI((k-6*res));
        B2(8,8) =  unitI((k-6*res));
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=8*res
        % keep u2,...,u8, fill u9
        B1(9,10) =  unitI((k-7*res));
        B2(9,9) =  unitI((k-7*res));
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    elseif k<=9*res
        % keep u2,...,u9, fill u10
        B1(10,9) = -unitI((k-8*res));
        B2(10,10) =  unitI((k-8*res));
        curv_rank_Stiefel_cplot(k)  = seccurv_Stiefel_canon(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eplot(k)  = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1);
        curv_rank_Stiefel_eTplot(k) = seccurv_Stiefel_euclid(zeros(n),B2,zeros(n),B1');
        curv_rank_Grassmann_plot(k) = seccurv_Grassmann(B2,B1);
        curv_rank_SOn_plot(k)       = ...
            seccurv_SOn([[zeros(n), -B2'];[B2, zeros(n)]],[[zeros(n), -B1'];[B1, zeros(n)]]);
    end
end

% display the results
figure;

plot(1:9*res, curv_rank_Grassmann_plot, 'b-.', ...
     1:9*res, curv_rank_Stiefel_cplot,  'k',...
     1:9*res, curv_rank_Stiefel_eplot,  'k:',...
     1:9*res, curv_rank_Stiefel_eTplot, 'k--',...
     1:9*res, curv_rank_SOn_plot, 'r-', ...
     ...
     LineWidth=2)
ax = gca;
ax.FontSize = 24; 
legend("Grassmann", "Stiefel (can.)", "Stiefel (Euclid)", "Stiefel (Euclid, $B_1^T$)",  "SO(n)", 'interpreter','latex')
title('Sectional curvature under increasing rank')
xlabel('$u$', 'interpreter','latex','FontSize',24)
xticks([0 res 2*res 3*res 4*res 5*res 6*res 7*res 8*res, 9*res])
xticklabels({'1', '2','3','4','5','6','7','8', '9', '10'})
ylabel('$K(0,B_1(u), 0, B_2(u))$', 'interpreter','latex','FontSize',24)
