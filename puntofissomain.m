%% P2-P0 FEM for Navier-Stokes
% GLOBAL DOFS ORDERING: first DoFs first component velocity (vertices and then edges);
%                       then DoFs second component (vertices and then edges)
%                       then pressure 
%% We use fix point iteration to cope with the nonlinearity
function [uh,pip] = puntofissomain(mesh,maxIter,nu,type);

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
nu = 1; 

%% SETTING THE QUADRATURE FORMULAS
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
A = sparse(2*(nver+nele),2*(nver+nele));
B = sparse(nver,2*(nver+nele));
b1 = zeros(nver+nele,1);   
b2 = zeros(nver+nele,1);  

vettarea = zeros(nver,1);  % needed later for the zero average constraint

%% we compute once and for all the matrices A and B, and the RHS
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
    vettarea([v1 v2 v3]) = vettarea([v1 v2 v3]) + (area/3)*ones(3,1);
    
    % setting local matrices
    KE = zeros(4,4);    
    BE = zeros(3,8);   

    %% quadrature formula
    for i=1:4
        for j=1:i-1
            KE(i,j) = KE(j,i);
        end
        for j=i:4
            for q=1:Nq
                tmp = JF*[xhq(q);yhq(q)]+[x1;y1];
                xq = tmp(1); yq = tmp(2);
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
    
    %% construct the mixed matrix

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
    
    % ordering of the edges
    %l1 = edges(iele,1); l2 = edges(iele,2); l3 = edges(iele,3);
    
     % array of the local2global DoFs
    dofg = [v1 v2 v3 nver+iele];          % 1st component
    dofgg = [dofg , dofg + (nver+nele)];  % add 2nd component
    
    %% update stiffness matrix
    A(dofgg,dofgg) = A(dofgg,dofgg) + KE;
    
    %% ASSEMBLE MIXED MATRIC
    dofp = [v1 v2 v3];    % GdL globali della pressione (indici dei vertici)
    B(dofp,dofgg) = B(dofp,dofgg) + BE;


%% quadrature for loading term
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
b = [b1;b2];

%% homogeneous BCs
internV = setdiff([1:1:nver+nele],boundary);
NL  = [internV,(nver+nele) + internV];
Ah = A(NL,NL);
Bh = B(:,NL);
fh = b(NL);
% clear A
% clear B
% clear b

%% Construct the global matrix and rhs
Kh = [Ah,Bh';Bh,zeros(nver)];
fh = [fh;zeros(nver,1)];


%% If the BCs are homogeneous Dirichlet, impose zero average constraint with multiplier
% OTHERWISE: comment the lines below
N = length(NL);
Kh = [Kh,[zeros(N,1);vettarea];zeros(1,N),vettarea',0];
fh = [fh;0];

%% NEW PART (fixed point iteration)
% The term stemming from the trilinear form is computed at each iteration

uh = zeros(2*(nver+nele),1);
max_iter = 1000;
n_iter = 1;
pip=1; %tolleranza

while (pip > 1.e-10 && n_iter <= maxIter) % arrest under two conditions
    n_iter;
    C = sparse(2*(nver+nele),2*(nver+nele));  % matrix C associated with the trilinear term
    
    uh1 = uh(1:length(uh)/2); uh2 = uh(length(uh)/2+1:end);
    
    for iele=1:nele
        CE = zeros(4,4);  % ``scalar'' block
    
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

        % local2global DoFs
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
                    
                    for ii=1:4
                        un1 = un1 + uT1(ii)*phihq(ii,q);
                        un2 = un2 + uT2(ii)*phihq(ii,q); 
                    end
                    tmp = phihq(i,q)*dot((JFIT*[gphihqx(j,q);gphihqy(j,q)]),[un1,un2]); 
                    CE(i,j) = CE(i,j) + tmp*whq(q); 
                end
                CE(i,j) = 2*area*CE(i,j);
            end
        end
        CE = [CE,zeros(4);zeros(4),CE];
        CE = 0.5*(CE - CE');
        C(dofgg,dofgg) = C(dofgg,dofgg) + CE;
    end

    %C = (C-C')/2;  % anti-simmetrize  | OPTIONAL AND COMMENTABLE |

    Ch = C(NL,NL);  % internal DoFs
    clear C

    %% update matrix with trilinear form
    Kh = [Ah+Ch,Bh';Bh,zeros(nver)];

    %% imposing zero average multiplier
    Kh = [Kh,[zeros(N,1);vettarea];zeros(1,N),vettarea',0];

    %% update step
    uhold=uh; 
    
    %% solve system and extract
    uh = zeros(2*(nver+nele),1);
    solh = Kh\fh;           
    uh(NL) = solh(1:N);            % estraggo le velocita'
    ph  = solh(N+1:N+nver);        % estraggo le pressioni

    %% update quantities to stop the iteration
    disp('Pip di punto fisso')
    pip = norm(uh-uhold)/norm(uh)
    disp('------------------------------')
    n_iter = n_iter + 1;
end

