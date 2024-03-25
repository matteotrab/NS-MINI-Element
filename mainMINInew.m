%% P2-P0 FEM for Navier-Stokes
% GLOBAL DOFS ORDERING: first DoFs first component velocity (vertices and then edges);
%                       then DoFs second component (vertices and then edges)
%                       then pressure 
%% We use fix point iteration to cope with the nonlinearity
function [errH1,errPR,errL2, diameter, n_iter] = mainMINInew(mesh, nu, TOL,plot,type);

iterBackTracking=0;

%%%%% Scelta forma da usare
isTilde=1; %uguale a 1 per usare la forma tildata, 0 altrimenti
isInverted=0; %uguale a 1 per usare la forma con la nuova definizione delle slide

carico = strcat('carico', num2str(type));
uexact = strcat('uexact', num2str(type));
uexactG = strcat('uexactG', num2str(type));
pexact = strcat('pexact', num2str(type));

%% extract geometric information
xv=mesh.xv;
yv=mesh.yv;
vertices=mesh.vertices;
edges=mesh.edges;
endpoints=mesh.endpoints;
boundary=mesh.boundary;
boundedges=mesh.boundedges;

%% viscosity
%nu = 1; 

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
% b = [b1;b2];

%% homogeneous BCs
internV = setdiff([1:1:nver+nele],boundary);
NL  = [internV,(nver+nele) + internV];
Ah = A(NL,NL);
Bh = B(:,NL);
% fh = b(NL);
% clear A
% clear B
% clear b

%% Construct the global matrix and rhs
Kh = [Ah,Bh';Bh,zeros(nver)];
% fh = [fh;zeros(nver,1)];


%% If the BCs are homogeneous Dirichlet, impose zero average constraint with multiplier
% OTHERWISE: comment the lines below
N = length(NL);
Kh = [Kh,[zeros(N,1);vettarea];zeros(1,N),vettarea',0];
% fh = [fh;0];

%% NEW PART (newton iteration)
% The term stemming from the trilinear form is computed at each iteration

uh = zeros(2*(nver+nele),1);

max_iter = 100;
n_iter = 1;
pip=1; %tconfronto iterate

if iterBackTracking>0
    [uh,pip] = puntofissomain(mesh,iterBackTracking, nu, type);
end

b1original=b1;
b2original=b2;

while (pip > TOL && n_iter < max_iter) % arrest under two conditions
    n_iter
    C = sparse(2*(nver+nele),2*(nver+nele));  % matrix C associated with the trilinear term
    
    uh1 = uh(1:length(uh)/2); 
    uh2 = uh(length(uh)/2+1:end);
    
    
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
        if isInverted==0
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
            if isTilde==1
                CE = 0.5*(CE - CE');
            end

       elseif isInverted==1
            %Pezzo uno a
            uno=zeros(8);
            unoa=zeros(8);
            unob=zeros(8);
            for i=1:8
                for j=1:8
                    for q=1:Nq
                        un1  = 0; un2  = 0;
                        %gradun1  = 0; gradun2  = 0;

                        for ii=1:4
                                un1 = un1 + uT1(ii)*phihq(ii,q);
                                un2 = un2 + uT2(ii)*phihq(ii,q); 


                        end

                    if i<4.5
                        if j<4.5
                            unoa(i,j)=unoa(i,j)+whq(q)*un1*dot((JFIT*[gphihqx(j,q);gphihqy(j,q)]),[1,0])*phihq(i,q);
                        elseif j>4.5
                            unoa(i,j)=unoa(i,j)+whq(q)*un2*dot((JFIT*[gphihqx(j-4,q);gphihqy(j-4,q)]),[1,0])*phihq(i,q);
                        end
                    elseif i>4.5
                        if j<4.5                                             
                            unoa(i,j)=unoa(i,j)+whq(q)*un1*dot((JFIT*[gphihqx(j,q);gphihqy(j,q)]),[0,1])*phihq(i-4,q);
                        elseif j>4.5                       
                            unoa(i,j)=unoa(i,j)+whq(q)*un2*dot((JFIT*[gphihqx(j-4,q);gphihqy(j-4,q)]),[0,1])*phihq(i-4,q);
                        end
                    end
                    end
                end
            end
            uno=2*area*unoa;
            if isTilde==1
                CE=(uno-uno')/2;
            elseif isTilde==0
                CE=uno;
            end
        
        end
        
        
        %Pezzo 2 a e b
        due=zeros(8);
        duea=zeros(8);
        dueb=zeros(8);
        for i=1:8
            for j=1:8
                for q=1:Nq
                    un1  = 0; un2  = 0;
                    %gradun1  = 0; gradun2  = 0;
%                     gradun1x =0;
%                     gradun2x =0;
%                     gradun1y =0;
%                     gradun2y =0;
                    gradun2 = [0;0];
                    gradun1 = [0;0];

                    for ii=1:4
                            un1 = un1 + uT1(ii)*phihq(ii,q);
                            un2 = un2 + uT2(ii)*phihq(ii,q); 
%                             gradun1x = gradun1x + uT1(ii)*gphihqx(ii,q);
%                             gradun2x = gradun2x + uT2(ii)*gphihqx(ii,q);
%                             gradun1y = gradun1y + uT1(ii)*gphihqy(ii,q);
%                             gradun2y = gradun2y + uT2(ii)*gphihqy(ii,q);

                            gradun1 = gradun1 + uT1(ii)*(JFIT*[gphihqx(ii,q);gphihqy(ii,q)]);
                            gradun2 = gradun2 + uT2(ii)*(JFIT*[gphihqx(ii,q);gphihqy(ii,q)]);
                    end
%                     gradun1=[gradun1x;gradun1y];
%                     gradun2=[gradun2x;gradun2y];
                    
                if i<4.5
                    if j<4.5
                        if isInverted==0
                            duea(i,j)=duea(i,j)+whq(q)*gradun1(1)*phihq(j,q)*phihq(i,q);
                            dueb(i,j)=dueb(i,j)+whq(q)*(JFIT(1,1:2)*...
                               [gphihqx(i,q);gphihqy(i,q)])*phihq(j,q)*un1;
                        elseif isInverted==1
                            duea(i,j)=duea(i,j)+whq(q)*gradun1(1)*phihq(j,q)*phihq(i,q);
                            dueb(i,j)=dueb(i,j)+whq(q)*dot((JFIT*[gphihqx(i,q);gphihqy(i,q)]),[un1,un2])*phihq(j,q);
                        end
                    elseif j>4.5
                        if isInverted==0
                            duea(i,j)=duea(i,j)+whq(q)*gradun1(2)*phihq(j-4,q)*phihq(i,q);
                            dueb(i,j)=dueb(i,j)+whq(q)*(JFIT(2,1:2)*...
                              [gphihqx(i,q);gphihqy(i,q)])*phihq(j-4,q)*un1;
                        elseif isInverted==1
                            duea(i,j)=duea(i,j)+whq(q)*gradun2(1)*phihq(j-4,q)*phihq(i,q);
                            dueb(i,j)=dueb(i,j)+0;
                        end
                    end
                elseif i>4.5
                    if j<4.5
                        if isInverted==0
                            duea(i,j)=duea(i,j)+whq(q)*gradun2(1)*phihq(j,q)*phihq(i-4,q);
                            dueb(i,j)=dueb(i,j)+whq(q)*(JFIT(1,1:2)*...
                              [gphihqx(i-4,q);gphihqy(i-4,q)])*phihq(j,q)*un2;
                        elseif isInverted==1
                            duea(i,j)=duea(i,j)+whq(q)*gradun1(2)*phihq(j,q)*phihq(i-4,q);
                            dueb(i,j)=dueb(i,j)+0;
                        end
                    elseif j>4.5 
                        if isInverted==0
                            duea(i,j)=duea(i,j)+whq(q)*gradun2(2)*phihq(j-4,q)*phihq(i-4,q);
                            dueb(i,j)=dueb(i,j)+whq(q)*(JFIT(2,1:2)*...
                              [gphihqx(i-4,q);gphihqy(i-4,q)])*phihq(j-4,q)*un2;
                        elseif isInverted==1
                            duea(i,j)=duea(i,j)+whq(q)*gradun2(2)*phihq(j-4,q)*phihq(i-4,q);
                            dueb(i,j)=dueb(i,j)+whq(q)*dot((JFIT*[gphihqx(i-4,q);gphihqy(i-4,q)]),[un1,un2])*phihq(j-4,q);
                        end
                    end
                end
                end
            end
        end
        
        %moltiplico per 1/2
        if isTilde==1
            due=area*(duea-dueb);
            %due=area*(duea-duea')/2;
        elseif isTilde==0
            due=2*area*duea;
        end
        
        trea=zeros(8,1);
        treb=zeros(8,1);
        tre=zeros(8,1);
        for i=1:8
            for q=1:Nq
                    un1  = 0; un2  = 0;
%                     gradun1x =0;
%                     gradun2x =0;
%                     gradun1y =0;
%                     gradun2y =0;
                    gradun2 = [0;0];
                    gradun1 = [0;0];
                    for ii=1:4
                            un1 = un1 + uT1(ii)*phihq(ii,q);
                            un2 = un2 + uT2(ii)*phihq(ii,q); 
%                             gradun1x = gradun1x + uT1(ii)*gphihqx(ii,q);
%                             gradun2x = gradun2x + uT2(ii)*gphihqx(ii,q);
%                             gradun1y = gradun1y + uT1(ii)*gphihqy(ii,q);
%                             gradun2y = gradun2y + uT2(ii)*gphihqy(ii,q);

                            gradun1 = gradun1 + uT1(ii)*(JFIT*[gphihqx(ii,q);gphihqy(ii,q)]);
                            gradun2 = gradun2 + uT2(ii)*(JFIT*[gphihqx(ii,q);gphihqy(ii,q)]);

                    end
%                     gradun1=[gradun1x;gradun1y];
%                     gradun2=[gradun2x;gradun2y];
                if i<4.5
                    if isInverted==0
                       trea(i,1)=trea(i,1)+whq(q)*dot(gradun1,[un1;un2])*...
                          phihq(i,q);
                       aaaaa=JFIT*[gphihqx(i,q);gphihqy(i,q)];
                       treb(i,1)=treb(i,1)+whq(q)*dot(aaaaa,[un1;un2])*...
                          un1;
                    elseif isInverted==1
                       trea(i,1)=trea(i,1)+whq(q)*(un1*gradun1x+un2*gradun2x)*phihq(i,q);
                       treb(i,1)=treb(i,1)+whq(q)*un1*dot((JFIT*[gphihqx(i,q);gphihqy(i,q)]),[un1,un2]);
                    end
                elseif i>4.5
                    if isInverted==0
                       trea(i,1)=trea(i,1)+whq(q)*dot(gradun2,[un1;un2])*...
                          phihq(i-4,q);
                       AAAAA=JFIT*[gphihqx(i-4,q);gphihqy(i-4,q)];
                       treb(i,1)=treb(i,1)+whq(q)*dot(AAAAA,[un1;un2])*...
                          un2;
                    elseif isInverted==1
                       trea(i,1)=trea(i,1)+whq(q)*(un1*gradun1y+un2*gradun2y)*phihq(i-4,q);
                       treb(i,1)=treb(i,1)+whq(q)*un2*dot((JFIT*[gphihqx(i-4,q);gphihqy(i-4,q)]),[un1,un2]);
                    end
                end
            end
        end
        
        if isTilde==1
            tre=area*(trea-treb);
            %tre=area*(trea-trea')/2;
        elseif isTilde==0
            tre=2*area*trea;
        end
                
        b1=b1original;
        b2=b2original;
        
        C(dofgg,dofgg) = C(dofgg,dofgg) + CE + due;
        b1(dofg) = b1(dofg) + tre(1:4,1);
        b2(dofg) = b2(dofg) + tre(5:8,1);
        
    
    end
    b=[b1;b2];
    fh = b(NL);
    fh = [fh;zeros(nver,1)];
    fh = [fh;0];

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
    Kh = sparse(Kh);
    solh = Kh\fh;           
    uh(NL) = solh(1:N);            % estraggo le velocita'
    ph  = solh(N+1:N+nver);        % estraggo le pressioni

    %% update quantities to stop the iteration
    pip = norm(uh-uhold)/norm(uh);
    n_iter = n_iter + 1;
end

if(n_iter == max_iter)
    disp('fixed point does not converge')
end    

uh1 = uh(1:length(uh)/2); 
uh2 = uh(length(uh)/2+1:end);

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

%% PLOT SOLUTIONS
if plot == 1
    name = 'Newton';
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