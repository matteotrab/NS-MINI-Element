clc
clear all
addpath ../../Meshes/  % serve per richiamare il file meshtrans.m
% in alternativa metterlo nella stessa cartella


H = 1./[2 4 8 16 32];

mesh1 = meshgen(H(1));
mesh2 = meshgen(H(2));
mesh3 = meshgen(H(3));
mesh4 = meshgen(H(4));
mesh5 = meshgen(H(5));


NDOFS = zeros(size(H));

NDOFS(1) = size(mesh1.vertices,1)+2*(size(mesh1.xv,1)+size(mesh1.endpoints,1)-length(mesh1.boundary)-length(mesh1.boundedges)); 
NDOFS(2) = size(mesh2.vertices,1)+2*(size(mesh2.xv,1)+size(mesh2.endpoints,1)-length(mesh2.boundary)-length(mesh2.boundedges));
NDOFS(3) = size(mesh3.vertices,1)+2*(size(mesh3.xv,1)+size(mesh3.endpoints,1)-length(mesh3.boundary)-length(mesh3.boundedges));
NDOFS(4) = size(mesh4.vertices,1)+2*(size(mesh4.xv,1)+size(mesh4.endpoints,1)-length(mesh4.boundary)-length(mesh4.boundedges));
NDOFS(5) = size(mesh5.vertices,1)+2*(size(mesh5.xv,1)+size(mesh5.endpoints,1)-length(mesh5.boundary)-length(mesh5.boundedges));



errH1 = zeros(size(H));
errPR = zeros(size(H));
errL2 = zeros(size(H));


[errH1(1),errPR(1),errL2(1)] = main_fix(mesh1);
[errH1(2),errPR(2),errL2(2)] = main_fix(mesh2);
[errH1(3),errPR(3),errL2(3)] = main_fix(mesh3);
[errH1(4),errPR(4),errL2(4)] = main_fix(mesh4);
[errH1(5),errPR(5),errL2(5)] = main_fix(mesh5);



figure(2)
loglog(H, errH1, 'LineStyle', '--', 'Marker','square', 'Linewidth', 2, 'MarkerSize', 12, 'Color', 'b')
hold on
loglog(H, errPR, 'LineStyle', '--', 'Marker','x', 'Linewidth', 2, 'MarkerSize', 12, 'Color', 'r')
loglog(H, errL2, 'LineStyle', '--', 'Marker','o', 'Linewidth', 2, 'MarkerSize', 8, 'Color', 'g')

loglog(H, H, 'LineStyle', '-', 'Linewidth', 2, 'Color', 'k')
loglog(H, H.^2, 'LineStyle', '-', 'Linewidth', 2, 'Color', 'k')

L = legend('$H^1$-\texttt{error velocity}', '$L^2$-\texttt{error pressure}', '$L^2$-\texttt{error velocity}', '$h$', '$h^2$');
set(L, 'Interpreter', 'latex', 'FontSize', 24, 'Location','southeast');
xlabel('$h$', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('errors', 'Interpreter', 'latex', 'FontSize', 20)



figure(3)
loglog(NDOFS, errH1, 'LineStyle', '--', 'Marker','square', 'Linewidth', 2, 'MarkerSize', 12, 'Color', 'b')
hold on
loglog(NDOFS, errPR, 'LineStyle', '--', 'Marker','x', 'Linewidth', 2, 'MarkerSize', 12, 'Color', 'r')
loglog(NDOFS, errL2, 'LineStyle', '--', 'Marker','o', 'Linewidth', 2, 'MarkerSize', 8, 'Color', 'g')
loglog(NDOFS, NDOFS.^(-0.5), 'LineStyle', '-', 'Marker','o', 'Linewidth', 2, 'Color', 'k')
loglog(NDOFS, NDOFS.^(-1), 'LineStyle', '-', 'Marker','o', 'Linewidth', 2, 'Color', 'k')

M = legend('$H^1$-\texttt{error velocity}','$L^2$-\texttt{error pressure}','$L^2$-\texttt{error velocity}', '$\texttt{NDOFD}^{-0.5}$', '$\texttt{NDOFD}^{-1}$');
set(M, 'Interpreter', 'latex', 'FontSize', 24, 'Location','northeast');
xlabel('\texttt{NDOFS}', 'Interpreter', 'latex', 'FontSize', 20)
ylabel('errors', 'Interpreter', 'latex', 'FontSize', 20)