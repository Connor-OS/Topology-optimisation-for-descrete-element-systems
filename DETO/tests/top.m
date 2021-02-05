% A DISCRETE ELEMENT TOPOLOGY OPTIMIZATION CODE BY CONNOR O'SHAUGHNESSY 2020%
nelx = 75;
nely = 25;
mass = 0.6;
Diam = 1;
rmin = 1.1;
kspr = 100;

function top(nelx,nely,mass,Diam,rmin,kspr)
tic
% INITIALIZEATION
m = mass*ones((nely*nelx)-floor(nely/2),1);
cut = 1.01*Diam;
xi = zeros(length(m),1);  yi = xi;
Fxe = zeros(length(m),1);   Fye = Fxe;
Fye(length(m)-floor(nelx/2)) = -0.1;
% Initial geometry
for j = 1:nely
    for i = 1:nelx-mod(j+1,2)
        xi((j-1)*nelx+i-floor((j-1)/2)) = Diam*(i-0.5)+(1-mod(j,2))*Diam/2-(nelx*Diam)/2;
        yi((j-1)*nelx+i-floor((j-1)/2)) = Diam*((j-1)*cos(pi/6)+0.0);
    end
end
% Neighbour list and equilibrium distances
neighbour_list
% Variables and plotting
loop = 0;
change = 1.;
x = xi; y = yi;
xy = sprintf('../dump/xy_%i_%i_%1.2f_%1.2f_%1.2f_%i.txt',nelx,nely,mass,Diam,rmin,kspr);
it = sprintf('../dump/it_%i_%i_%1.2f_%1.2f_%1.2f_%i.txt',nelx,nely,mass,Diam,rmin,kspr);
XY = fopen(xy,'w');
Itt = fopen(it,'w');
% START ITERATION
while change > 0.004
    loop = loop+1;
    mold = m;
    for i = 1:length(m)
        for s = 1:nn(i)
            j = N(i,s);
            k(i,s) = m(i)^2 *m(j)^2  * kspr;
        end
    end
    % DEM-ANALYSIS assuming harmonic potential
    Qmin
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
        Mesh_filter
    end
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    OC
    % PER PARTICLE STRESS
    stress_perpart
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
post_processing
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
