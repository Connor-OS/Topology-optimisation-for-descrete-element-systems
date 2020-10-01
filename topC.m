% A DISCRETE ELEMENT TOPOLOGY OPTIMIZATION CODE BY CONNOR O'SHAUGHNESSY 2020%
%%% Specific function call here for rocket application

function topC(nelx,nely,mass,Diam,rmin,kspr)
tic
% INITIALIZEATION
m = mass*ones((nely*nelx)-floor(nely/2),1);
cut = 1.01*Diam;
xi = zeros(length(m),1);  yi = xi;
Fxe = zeros(length(m),1);   Fye = Fxe;
Fye(length(m)-floor(nelx/2)) = -1;
% Initial geometry
for j = 1:nely
    for i = 1:nelx-mod(j+1,2)
        xi((j-1)*nelx+i-floor((j-1)/2)) = Diam*(i-0.5)+(1-mod(j,2))*Diam/2-(nelx*Diam)/2;
        yi((j-1)*nelx+i-floor((j-1)/2)) = Diam*((j-1)*cos(pi/6)+0.0);
    end
end
% Neighbour list and equilibrium distances
[N,Nf,nn,nnf,Li,Lif] = nlist(xi,yi,cut,rmin,Diam);
% Variables and plotting
loop = 0;
change = 1.;
x = xi; y = yi;
xy = sprintf('xy_%i_%i_%1.2f_%1.2f_%1.2f_%i.txt',nelx,nely,mass,Diam,rmin,kspr);
it = sprintf('it_%i_%i_%1.2f_%1.2f_%1.2f_%i.txt',nelx,nely,mass,Diam,rmin,kspr);
XY = fopen(xy,'w');
Itt = fopen(it,'w');
% START ITERATION
while change > 0.002
    loop = loop+1;
    mold = m;
    for i = 1:length(m)
        for s = 1:nn(i)
            j = N(i,s);
            k(i,s) = m(i)^2 *m(j)^2  * kspr;
        end
    end
    % DEM-ANALYSIS assuming harmonic potential
    [x,y,L] = DEM(x,y,nn,N,nelx,nely,k,Li,Fxe,Fye,Diam,kspr,loop);
    % Cost function and sensitivities
    Eer = 0;
    dc = zeros(length(xi),1);
    for i = 1:length(xi)
        for s = 1:nn(i)
            j = N(i,s);
            Eer = Eer+1/4*m(i)^2*m(j)^2*kspr*(L(i,s)-Li(i,s))^2; % 1/4 cause double counting each interaction
            dc(i) = dc(i)-0.5*m(i)*m(j)^2*kspr*(L(i,s)-Li(i,s))^2;
        end
    end
    % Filtering of sensitivities (aka coarse graining)
    if rmin > Diam
        [dc] = check(rmin,m,Nf,nnf,dc,Lif);
    end
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [m] = OC(m,mass,dc);
    % PER PARTICLE STRESS
    [W] = strs(m,N,nn,k,L,Li,x,y,Diam);
    % PRINT RESULTS
    change = max(abs(m-mold));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%6.4f',Eer) ...
       ' Mass.: ' sprintf('%6.3f ',sum(m)/(length(xi))) ...
       ' ch.:' sprintf('%6.3f',change) ' Time.: ' sprintf('%4.2f',toc)])
    % PLOT DENSITIES
    cc=[1-m,1-m,1-m];
    linkdata on
    scatter(x,y,Diam*50,cc,'filled')
    % WRITE DUMP FILES
    fprintf(XY,[sprintf('%4i\n',length(m)),'Frame.: ',sprintf('%4i\n',loop)]);
    for i=1:length(x)
        fprintf(XY,[sprintf('%4i ',i) sprintf('%6.3f ',x(i)) sprintf('%6.3f ',y(i)) ...
            sprintf('%6.3f ',m(i)) sprintf('%6.3f ',W(i,1)) sprintf('%6.3f ',W(i,2)) ...
            sprintf('%6.3f ',W(i,3)) sprintf('%6.3f ',W(i,4)) sprintf('%6.3f ',W(i,5))...
            sprintf('%6.3f\n',W(i,6))]);
    end
    fprintf(Itt,[' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',Eer) ...
        ' Mass.: ' sprintf('%6.3f',sum(m)/(length(xi))) ...
        ' ch.: ' sprintf('%6.3f',change ) ' Time.: ' sprintf('%6.3f\n',toc)]);
 end
 fprintf(Itt,[' Sim time.: ' sprintf( '%6.3f\n', toc)]);
%%%%% POST PROCESSING %%%%%%%
p1 = 0; p2 = 1;
while p2-p1 > 0.0001
    p = (p1+p2)/2;
    mfinal = ceil(m-p);
    if sum(mfinal)-mass*length(m) > 0
        p1 = p;
    else
        p2 = p;
    end
end
m = mfinal;
for i = 1:length(m)
    for s = 1:nn(i)
        j = N(i,s);
        k(i,s) = m(i)^2 *m(j)^2  * kspr;
    end
end
cc = [1-m,1-m,1-m];
scatter(x, y, Diam*50, cc, 'filled')
[x,y,L] = DEM(x,y,nn,N,nelx,nely,k,Li,Fxe,Fye,Diam,kspr,loop);
[W] = strs(m,N,nn,k,L,Li,x,y,Diam);
toc
for i = 1:length(xi)
    for s = 1:nn(i)
        Eer = Eer+1/4*m(i)^2*m(j)^2*kspr*(L(i,s)-Li(i,s))^2;
    end
end
fprintf(XY,[sprintf('%4i\n',length(m)),'Frame.: ',sprintf('%4i\n',loop)]);
for i=1:length(x)
    fprintf(XY,[sprintf('%4i ',i) sprintf('%6.3f ',x(i)) sprintf('%6.3f ',y(i)) ...
        sprintf('%6.3f ',m(i)) sprintf('%6.3f ',W(i,1)) sprintf('%6.3f ',W(i,2)) ...
        sprintf('%6.3f ',W(i,3)) sprintf('%6.3f ',W(i,4)) sprintf('%6.3f ',W(i,5))...
        sprintf('%6.3f\n',W(i,6))]);
end
fprintf(Itt,[' Final Obj.: ' sprintf( '%10.4f', Eer)]);
fclose('all');
%%%%%%%%%% Neighbour lists and initial lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,Nf,nn,nnf,Li,Lif]=nlist(xi,yi,cut,rmin,Diam)
nn = zeros(length(xi),1);
nnf = zeros(length(xi),1);
N = zeros(length(xi),6); Li=N;
for i = 1:length(xi)
    for j = 1:length(xi)
        L = sqrt( (yi(i)-yi(j))^2 + (xi(i)-xi(j))^2 );
        if i ~= j && L < cut
            nn(i) = nn(i)+1;
            N(i,nn(i)) = j; Li(i,nn(i)) = L;
        end
        if i~=j && L < rmin
            nnf(i) = nnf(i)+1;
            Nf(i,nnf(i)) = j; Lif(i,nnf(i)) = L;
        end
    end
end
if rmin < Diam
    Nf = 0;
    Lif = 0;
end
%%%%%%%%%% DEM energy minimisation by steepest descent %%%%%%%%%%%%%%%%%%%%
function [x,y,L] = DEM(x,y,nn,N,nelx,nely,k,Li,Fxe,Fye,Diam,kspr,loop)
damp = max(300,3*kspr);
tol = Diam;
L=Li;
while tol > Diam*nely/(kspr*2e7)
    Fx=Fxe;  Fy = Fye;
    for i=1:length(x)
        for s=1:nn(i)
            j = N(i,s);
            dx = x(j)-x(i);     dy = y(j)-y(i);
            L(i,s) = sqrt( dx^2 + dy^2 );
            if j > i
                F = k(i,s) * (L(i,s) - Li(i,s));
                Fxx = F * dx/L(i,s);    Fyy = F * dy/L(i,s);
                Fx(i) = Fx(i) + Fxx;    Fx(j) = Fx(j) - Fxx;
                Fy(i) = Fy(i) + Fyy;    Fy(j) = Fy(j) - Fyy;
            end
        end
    end
    % constraints and displacement
    Fy(1) = 0;
    Fy(nelx) = 0;
    x = x + Fx/damp;
    y = y + Fy/damp;
    % criterion
    tol = max(abs(Fx)+abs(Fy))/damp;
end
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mnew] = OC(m,mass,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-8)
    lmid = 0.5*(l2+l1);
    mnew = max(0,max(m-move,min(1.,min(m+move,m.*sqrt(-dc./lmid)))));
    if sum(mnew) - mass * length(m) > 0
        l1 = lmid;
    else
        l2 = lmid;
    end
end
%%%%%%%%%% PER PARTICLE STRESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [W] = strs(m,N,nn,k,L,Li,x,y,Diam)
W = zeros(length(m),6);
for i = 1:length(m)
    for s = 1:nn(i)
        j = N(i,s);
        F = k(i,s)*(L(i,s)-Li(i,s));
        dx = x(j) - x(i);       dy = y(j) - y(i);
        Fx1 = F * dx/L(i,s);    Fy1 = F * dy/L(i,s);
        Fx2 = -Fx1;             Fy2 = -Fy1;
        W(i,1) = W(i,1)+0.5*(Fx1*x(i)+Fx2*x(j));
        W(i,2) = W(i,2)+0.5*(Fy1*y(i)+Fy2*y(j));
        W(i,3) = W(i,3)+0.5*(Fy1*x(i)+Fy2*x(j));
        W(i,4) = W(i,4)+0.5*(Fx1*y(i)+Fx2*y(j));
        W(i,5) = (W(i,1)+W(i,2))/2; % Hydrostatic stress
        W(i,6) = sqrt(0.5*((W(i,1)-W(i,2))^2+W(i,2)^2+W(i,1)^2)+3*W(i,3)^2); % Deviatoric
    end
end
pvol =((pi*Diam^2)/4)/0.9069;
W = (W/pvol);
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn] = check(rmin,m,Nf,nnf,dc,Lif)
dcn = (dc.*m.*m)*rmin;
for i = 1:length(m)
    tot = 0;
    for s = 1:nnf(i)
        j = Nf(i,s);
        fac = rmin-Lif(i,s);
        tot = tot+fac;
        dcn(i) = dcn(i)+fac*m(j)*m(j)*dc(j);
    end
    dcn(i) = dcn(i)/(m(i)*m(i)*tot);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Enrico Masoero and Connor O'Shaughnessy  %
% School of Engineering Newcastle University                               %
% Please sent your comments to the author: c.o'shaughnessy1@ncl.ac.uk      %
%                                                                          %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but does not guaranty that the code is   %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
