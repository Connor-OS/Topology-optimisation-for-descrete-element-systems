% A DISCRETE ELEMENT TOPOLOGY OPTIMIZATION CODE BY CONNOR O'SHAUGHNESSY 2020%
%%% Specific function call here for rocket application
function topDEM(nelx,nely,mass,Diam,rmin,kspr)
% INITIALIZE
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
loop = 0;
change = 1.;
x = xi; y = yi;
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
    [x,y,L] = DEM(x,y,nn,N,nelx,nely,k,Li,Fxe,Fye,Diam,kspr);
    % Cost function and sensitivities
    Eer = 0;
    dc = zeros(length(xi),1);
    for i = 1:length(xi)
        for s = 1:nn(i)
            j = N(i,s);
            Eer = Eer+1/4*m(i)^2*m(j)^2*kspr*(L(i,s)-Li(i,s))^2;
            dc(i) = dc(i)-m(j)^2*m(i)*kspr*(L(i,s)-Li(i,s))^2;
        end
    end
    % Filtering of sensitivities
    if rmin > Diam
        [dc] = check(rmin,m,Nf,nnf,dc,Lif);
    end
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [m] = OC(m,mass,dc);
    % PRINT RESULTS
    change = max(abs(m-mold));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%6.4f',Eer) ...
       ' Mass.: ' sprintf('%6.3f ',sum(m)/(length(xi))) ...
       ' ch.:' sprintf('%6.3f',change) ' Time.: ' sprintf('%4.2f',toc)])
    % PLOT DESIGN
    cc=[1-m,1-m,1-m];
    linkdata on
    scatter(x,y,Diam*50,cc,'filled')
end
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
    Nf = 0; Lif = 0;
end
%%%%%%%%%% DEM energy minimisation by steepest descent %%%%%%%%%%%%%%%%%%%%
function [x,y,L] = DEM(x,y,nn,N,nelx,nely,k,Li,Fxe,Fye,Diam,kspr)
damp = max(300,3*kspr);
tol = Diam;
L=Li;
while tol > Diam*nely/(kspr*2e6)
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
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn] = check(rmin,m,Nf,nnf,dc,Lif)
dcn = (dc.*m.*m)*rmin;
for i = 1:length(m)
    tot = 0;
    for s = 1:nnf(i)
        j = Nf(i,s);
        fac = rmin-Lif(i,s);
        tot = tot+fac;
        dcn(i) = dcn(i) + fac*m(j)*m(j)*dc(j);
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
