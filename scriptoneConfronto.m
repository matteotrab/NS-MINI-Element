
% SCRIPT PER NAVIER STOKES CON MINI, confronto metodi iterativi 

 clear all, close all
 load quadmeshes.mat

 %% SCELTA MESH
 meshes = [mesh;mesh1; mesh2; mesh3];
 endmesh = length(meshes)-1;

%% SCELTA VISCOSITA' E TOLLERANZA
nu = 1.5;
TOL = 10^(-10);

%% SCELTA DATI DEL PROBLEMA
% '1' per dati laboratorio, '2' per problema alternativo
type = 1;

carico = strcat('carico', num2str(type));
uexact = strcat('uexact', num2str(type));
uexactG = strcat('uexactG', num2str(type));
pexact = strcat('pexact', num2str(type));


%% PLOT CARICO, VELOCITA' E PRESSIONE ESATTE
figure()
[X,Y] = meshgrid(0:0.05:1,0:0.05:1);
surf(X,Y,feval(uexact,X,Y,1));
colorbar
title('Prima componente velocita esatta')


figure()
surf(X,Y,feval(uexact,X,Y,2))
title('Seconda componente velocita esatta')
colorbar

figure()
surf(X,Y,feval(pexact,X,Y))
title('Pressione esatta')
colorbar
view([12 60])

figure()
[X,Y] = meshgrid(0:0.05:1,0:0.05:1);
quiver(X,Y,feval(carico,X,Y,1,nu),feval(carico,X,Y,2,nu));
title('Carico quiver')

% autoArrangeFigures(2,2,1)
% drawnow

%% CICLO SULLE MESH

for j=1:endmesh
    j
    plot = 0;
    % per eventuale partenza da stokes
    %[uh1, uh2] = mainStokesMini(meshes2(j),nu);
    %uh = [uh1; uh2];

    % per plot solo di prima e ultima approssimazione
    if j == 1 || j == endmesh
        plot = 1;
    end
    [errH1IT, errPRIT, errL2IT, diameterIT, n_iterIT]=mainMINIfix(meshes(j),nu,TOL,plot,type);
    [errH1NEW, errPRNEW, errL2NEW, diameterNEW, n_iterNEW]=mainMINInew(meshes(j),nu, TOL,plot,type);


    h(j)         = diameterIT;

    errH1_vectIT(j)= errH1IT;
    errPR_vectIT(j)= errPRIT;
    nit_vectIT(j)  = n_iterIT;

    errH1_vectNEW(j)= errH1NEW;
    errPR_vectNEW(j)= errPRNEW;
    nit_vectNEW(j)  = n_iterNEW; 
end

%% PLOT ERRORI
% confronto con errore lineare sia nelle velocità che nelle pressioni
figure()
loglog(h,errH1_vectIT,'k');
title('Errore H^1 velocità')
subtitle('Iterata di punto fisso')
xlabel('h');
ylabel('H^1 error velocity')
hold on
loglog(8*[1e-1 1e-2],64*[1e-2 1e-3],'b--')
xlim([h(end), h(1)])
legend('errore','andamento lineare','Location','northwest')

figure()
loglog(h,errPR_vectIT,'k');
title('Errore L^2 pressioni')
subtitle('Iterata di Punto Fisso')
xlabel('h');
ylabel('L^2 error')
hold on
loglog(4*[1e-1 1e-2],16*[1e-2 1e-3],'g--')
xlim([h(end), h(1)])
legend('errore','andamento lineare','Location','northwest')

figure()
loglog(h,errH1_vectNEW,'k');
title('H^1 error velocity')
subtitle('Iterata di Newton')
xlabel('h');
ylabel('H^1 error velocity')
hold on
loglog(8*[1e-1 1e-2],64*[1e-2 1e-3],'b--')
xlim([h(end), h(1)])
legend('errore','andamento lineare','Location','northwest')


figure()
loglog(h,errPR_vectNEW,'k');
title('Errore L^2 pressioni')
subtitle('Iterata di Newton')
xlabel('h');
ylabel('L^2 error')
hold on
loglog(4*[1e-1 1e-2],16*[1e-2 1e-3],'g--')
xlim([h(end), h(1)])
legend('errore','andamento lineare','Location','northwest')


%% CONFRONTO ERRORI E NUMERO ITERATE

disp('Diametri successivi delle mesh: ')
h

disp('risultati iterata PUNTO FISSO al variare della mesh: ')

disp("Numero iterate: ")
nit_vectIT
disp('Errore velocità H1')
errH1_vectIT
disp('Errore pressioni L2')
errPR_vectIT

disp('risultati iterata NEWTON al variare della mesh: ')

disp("Numero iterate: ")
nit_vectNEW
disp('Errore velocità H1')
errH1_vectNEW
disp('Errore pressioni L2')
errPR_vectNEW

autoArrangeFigures();