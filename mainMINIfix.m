%% MINI FEM for Navier-Stokes con iterazione PUNTO FISSO
% GLOBAL DOFS ORDERING: first DoFs first component velocity (vertices and then edges);
%                       then DoFs second component (vertices and then edges)
%                       then pressure 
%% We use fix point iteration to cope with the nonlinearity
function [errH1,errPR,errL2, diameter, n_iter] = mainMINIfix(mesh, nu, TOL,plot,type);

%% extract geometric information
xv=mesh.xv;
yv=mesh.yv;
vertices=mesh.vertices;
edges=mesh.edges;
endpoints=mesh.endpoints;
boundary=mesh.boundary;
boundedges=mesh.boundedges;

carico = strcat('carico', num2str(type));
uexact = strcat('uexact', num2str(type));
uexactG = strcat('uexactG', num2str(type));
pexact = strcat('pexact', num2str(type));

%% viscosity
% nu = 1; 

%% SETTING THE QUADRATURE FORMULAS
% setto quadratura di ordine 4 perchè wc-scenario ho gradbolla.gradbolla
fdq = 'degree=7';
[xhq,yhq,whq]=quadratura(fdq);
Nq = length(xhq);

%% some measures
nver  = length(xv);         
nedge = size(endpoints,1);  
nele = size(vertices,1);    

%% basis functions at the quadrature nodes on the reference element
phihq = zeros(4,Nq);
for i=1:4
    for q=1:Nq
        phihq(i,q) = phih2(i,xhq(q),yhq(q));
    end
end

%% gradient of basis functions at the quadrature nodes on the reference element
gphihqx = zeros(4,Nq);
gphihqy = zeros(4,Nq);
for i=1:4
    for q=1:Nq
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end
end

%% initialize some parts of the global matrix
A = sparse(2*(nver+nele),2*(nver+nele));  % stiffness matrix
B = sparse(nver,2*(nver+nele));           % mixed matrix
b1 = zeros(nver+nele,1);                  % loading part 1
b2 = zeros(nver+nele,1);                  % loading part 2 


vettarea = zeros(nver,1);  % needed later for the zero average constraint


%% we compute once and for all the matrices A and B, and the RHS
% matrici che non necessitano di essere aggiornate durante iterate di punto
% fisso
for iele=1:nele
    %% recover element information

    % vertices number
    v1 = vertices(iele,1); v2 = vertices(iele,2); v3 = vertices(iele,3);
    
    % vertices coordinates
    x1 = xv(v1); y1 = yv(v1); x2 = xv(v2); y2 = yv(v2); x3 = xv(v3); y3 = yv(v3);
    
    % Jacobian
    JF = [x2-x1 x3-x1
        y2-y1 y3-y1];
    
    % inverse Jacobian
    JFI = inv(JF);
    
    % transpose inverse Jacobian
    JFIT = JFI';
    
    % areas
    area = 1/2*det(JF);
    %% CAPIRE PERCHE, PRESO DA MINI, IN STOKES USA COMMENTATO
    vettarea([v1 v2 v3]) = vettarea([v1 v2 v3]) + (area/3)*ones(3,1);
    %VettArea(iele)=area; % we need this later

    % setting local matrices
    % CONTROLLARE DIMENSIONI
    KE = zeros(4,4);    
    BE = zeros(3,8);   

    %% quadrature formula
    % pezzo forma bilineare grad grad
    for i=1:4
        for j=1:i-1
            KE(i,j) = KE(j,i);
        end
        for j=i:4
            for q=1:Nq
                tmp = JF*[xhq(q);yhq(q)]+[x1;y1];
                xq = tmp(1); 
                yq = tmp(2);
                tmp = dot(...
                      (JFIT*[gphihqx(j,q);gphihqy(j,q)]),...
                      (JFIT*[gphihqx(i,q);gphihqy(i,q)]));
                KE(i,j) = KE(i,j) + tmp*whq(q);
            end
            KE(i,j) = 2*area*KE(i,j);           
        end
    end
    
    %% duplicate scalar stiffness matrix
    KE = nu*[KE,zeros(4);zeros(4),KE];
    
    %% quadrature (mixed part)
    % la b coinvolge sia pressioni che velocità. Sfruttiamo le scalari per
    % costruirla appropiatamente
    for i=1:8  % 8 is the dimension of Vh
        for j=1:3  % 3 id the dimension of Qh
            for q=1:Nq
                % REMARK: the DIV of the first four basis functions for the velocities
                % ... is the x derivative of the scalar basis functions
                %         the DIV of the last four basis functions for the velocities
                % ... is the y derivative of the scalar basis functions
                if i<4.5 % first four basis elements
                    tmp = (JFIT(1,1:2)*[gphihqx(i,q);gphihqy(i,q)])*...
                        phihq(j,q);
                elseif i > 4.5 % last four basis elements
                    tmp = (JFIT(2,1:2)*[gphihqx(i-4,q);gphihqy(i-4,q)])*...
                        phihq(j,q);
                end
                BE(j,i)= BE(j,i) + tmp*whq(q);
            end
            BE(j,i) = 2*area*BE(j,i);
        end
    end
    
    %% ASSEMBLING THE GLOBAL MATRICES
        
    % array of the local2global DoFs
    dofg = [v1 v2 v3 nver+iele];          % 1st component
    dofgg = [dofg , dofg + (nver+nele)];  % add 2nd component

    %% update stiffness matrix
    A(dofgg,dofgg) = A(dofgg,dofgg) + KE;
    
    %% update mixed matrix
    dofp = [v1 v2 v3];    % GdL globali della pressione (indici dei vertici)
    B(dofp,dofgg) = B(dofp,dofgg) + BE;

    %% Compute the loading term, RHS
    % in navier stoker chiama function, confrontare
    fE1 = zeros(4,1);
    fE2 = zeros(4,1);

    for i=1:4
        for q=1:Nq
            tmp = JF*[xhq(q);yhq(q)]+[x1;y1];
            xq = tmp(1);
            yq = tmp(2);
            fE1(i) = fE1(i) + feval(carico,xq,yq,1,nu)*phihq(i,q)*whq(q);
            fE2(i) = fE2(i) + feval(carico,xq,yq,2,nu)*phihq(i,q)*whq(q);
        end
        fE1(i) = 2*area*fE1(i);
        fE2(i) = 2*area*fE2(i);
    end
    b1(dofg) = b1(dofg) + fE1;
    b2(dofg) = b2(dofg) + fE2;

end

b = [b1;b2]; %RHS?

% VEDERE FUNCTION E CONFRONTARE CON MINI MAIN
% % questo è quello che c'e nello stokes
% b1 = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints,1);
% b2 = P2load(xv,yv,vertices,edges,boundary,boundedges,endpoints,2);
% fh = [b1; b2];

%% homogeneous BCs
internV = setdiff([1:1:nver+nele],boundary);
NL  = [internV,(nver+nele) + internV];
Ah = A(NL,NL);
Bh = B(:,NL);
fh = b(NL);
clear A
clear B
clear b


%% Construct the global matrix and rhs
Kh = [Ah,Bh';Bh,zeros(nver)];
fh = [fh;zeros(nver,1)];


%% include the multiplier conditions
N = length(NL);
Kh = [Kh,[zeros(N,1);vettarea];zeros(1,N),vettarea',0];
fh = [fh;0];


%% NEW PART (fixed point iteration)
% The term stemming from the trilinear form is computed at each iteration

%%
 uh = zeros(2*(nver+nele),1); %partendo da 0
% uh = uhS; % partendo da Stokes

max_iter = 1000 ;
n_iter = 1;
%TOL = 10^(-10); % impostare tolleranza per criterio di arresto
conf_iter=1;


while (conf_iter > TOL && n_iter < max_iter) % arrest under two conditions
    n_iter;
    C = sparse(2*(nver+nele),2*(nver+nele));  % matrix C associated with the trilinear term
    % stessa dimensione di A 

    uh1 = uh(1:length(uh)/2); uh2 = uh(length(uh)/2+1:end);
    
    for iele=1:nele
        CE = zeros(4,4);  % blocco scalare 4x4
    
        %% recover element information

        % vertices number
        v1 = vertices(iele,1); v2 = vertices(iele,2); v3 = vertices(iele,3);

        % vertices coordinates
        x1 = xv(v1); y1 = yv(v1); x2 = xv(v2); y2 = yv(v2); x3 = xv(v3); y3 = yv(v3);

        % Jacobian
        JF = [x2-x1 x3-x1
            y2-y1 y3-y1];

        % inverse Jacobian
        JFI = inv(JF);

        % transpose inverse Jacobian
        JFIT = JFI';

        % areas
        area = 1/2*det(JF);

        % array of the local2global DoFs
        dofg = [v1 v2 v3 nver+iele];          % 1st component
        dofgg = [dofg , dofg + (nver+nele)];  % add 2nd component

        uT1 = uh1(dofg);  
        uT2 = uh2(dofg); 
        
        %% quadrature

        for i=1:4
            for j=1:4
                for q=1:Nq
                    % uh al current step
                    un1  = 0; un2  = 0;
                    
                    % u come combinazione lineare dei coefficienti ottenuti
                    % per le funzioni di base
                    for ii=1:4
                        un1 = un1 + uT1(ii)*phihq(ii,q);
                        un2 = un2 + uT2(ii)*phihq(ii,q); 
                    end
                    
                    % gradphi
                    tmp = phihq(i,q)*dot((JFIT*[gphihqx(j,q);gphihqy(j,q)]),[un1,un2]); 
                    CE(i,j) = CE(i,j) + tmp*whq(q); 
                end
                CE(i,j) = 2*area*CE(i,j);
            end
        end
        CE = [CE,zeros(4);zeros(4),CE]; % blocco scalare
        CE = 0.5*(CE - CE'); %costruisco a1 tildata antisimmetrica

        %% Aggiorno matrice associata a a tilde
        C(dofgg,dofgg) = C(dofgg,dofgg) + CE;
    end

    Ch = C(NL,NL);  % internal DoFs
    clear C

    %% update matrix with trilinear form
    Kh = [Ah+Ch,Bh';Bh,zeros(nver)];

    %% imposing zero average multiplier
    Kh = [Kh,[zeros(N,1);vettarea];zeros(1,N),vettarea',0];  

    %% update step
    uhold=uh; 
    
    %% solve system and extract velocity and pressure
    % dimensioni coerenti con MINI
    uh = zeros(2*(nver+nele),1);
    Kh = sparse(Kh);
    solh = Kh\fh;    


    uh(NL) = solh(1:N);            % estraggo le velocita'
    uh1 = uh(1:length(uh)/2);      % estraggo la prima componente
    uh2 = uh(length(uh)/2+1:end);  % estraggo la seconda componente
    ph  = solh(N+1:N+nver);        % estraggo le pressioni

    % solo per confronto, togliere

    %% update quantities to stop the iteration
    % criterio di arresto
    conf_iter = norm(uh-uhold)/norm(uh);
    n_iter = n_iter + 1;
end

if(n_iter == max_iter)
    disp('fixed point does not converge')
end    


%% QUADRATURE FORMULAS
fdq = 'degree=4';
[xhq,yhq,whq]=quadratura(fdq);
Nq = length(xhq);


%% basis functions at the quadrature nodes on the reference element
phihq = zeros(4,Nq);
for i=1:4
    for q=1:Nq
        phihq(i,q) = phih2(i,xhq(q),yhq(q));
        [gx,gy] = gradhphih2(i,xhq(q),yhq(q));
        gphihqx(i,q) = gx;
        gphihqy(i,q) = gy;
    end
end

%% Computing the H1 and L2 errors for the velocity, L2 error for the pressures
errL2sq = 0;   % errore L2 (velocit?)
errH1sq = 0;   % errore H1 (velocit?)
preL2sq = 0;   % errore L2 (pressioni) 

for iele=1:nele
    %% recover element information

    % vertices numbers
    v1 = vertices(iele,1); v2 = vertices(iele,2); v3 = vertices(iele,3);
    
    % vertices coordinates
    x1 = xv(v1); y1 = yv(v1); x2 = xv(v2); y2 = yv(v2); x3 = xv(v3); y3 = yv(v3);
    
    % Jacobian
    JF = [x2-x1 x3-x1
          y2-y1 y3-y1];

    % inverse Jacobian
    JFI = inv(JF);
    
    % transpose inverse Jacobian
    JFIT = JFI';

    % area
    area = 1/2*det(JF);

    %% array of the local2global DoFs (velocity and pressure)
    dofg = [v1 v2 v3 nver+iele];
    uT1 = uh1(dofg);     
    uT2 = uh2(dofg);     
    pT = ph([v1 v2 v3]);

    sq  = 0;
    sqH = 0;
    sqP = 0;
    
    %% quadrature
    for q=1:Nq
        tmp1  = 0;       % values of uh (first component)
        tmp2  = 0;       % values of uh (second component)
        tmpH1 = [0;0];   % values of gradient of uh (first component)
        tmpH2 = [0;0];   % values of gradient of uh (second component)
        tmpP = 0;        % valori della pressione ph

        for i=1:4
            tmp1 = tmp1 + uT1(i)*phihq(i,q);
            tmp2 = tmp2 + uT2(i)*phihq(i,q);
            tmpH1= tmpH1 + uT1(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]);
            tmpH2= tmpH2 + uT2(i)*(JFIT*[gphihqx(i,q);gphihqy(i,q)]);
        end
        tmpP = pT(1)*phihq(1,q) + pT(2)*phihq(2,q) + pT(3)*phihq(3,q);
        pos = JF*[xhq(q);yhq(q)]+[x1;y1];
        xq = pos(1);
        yq = pos(2);
        sq = sq + (feval(uexact,xq,yq,1)-tmp1)^2*whq(q) + (feval(uexact,xq,yq,2)-tmp2)^2*whq(q);
        sqH = sqH + norm(feval(uexactG,xq,yq,1) - tmpH1)^2*whq(q) ...
            + norm(feval(uexactG,xq,yq,2) - tmpH2)^2*whq(q);
        sqP = sqP + (feval(pexact,xq,yq)-tmpP)^2*whq(q);
    end
    %
    sq = sq*2*area;
    sqH = sqH*2*area;
    sqP = sqP*2*area;
    %
    errL2sq = errL2sq + sq;
    errH1sq = errH1sq + sqH;
    preL2sq = preL2sq + sqP;
    %
end

%% Errors
errL2 = sqrt(errL2sq); 
errH1 = sqrt(errH1sq);
errPR = sqrt(preL2sq); 

%% PLOT SOLUTIONS per prima e ultima mesh
if plot == 1
    name = 'Punto Fisso';
    plots_pressione_e_quiver
% pause()
end

%% SAVE DIAMETER OF THE MESH
diameter=0;
for j=1:length(mesh.endpoints)
    x1=mesh.xv(mesh.endpoints(j,1));
    y1=mesh.yv(mesh.endpoints(j,1));
    x2=mesh.xv(mesh.endpoints(j,2));
    y2=mesh.yv(mesh.endpoints(j,2));
    length_edge=pdist([x1 y1; x2 y2],'euclidean');
    if length_edge>diameter
        diameter=length_edge;
    end
end


