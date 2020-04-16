%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function topC(nelx,nely,mass,penal,rmin,Diam,cut,kspr,damp);
tic
% INITIALIZE
m = mass*ones((nely*nelx)-floor(nely/2),1);
xi = zeros(length(m),1);  yi=xi;
Fxe = zeros(length(m),1);   Fye=Fxe; 
Fye( length(m) -floor(nelx/2) ) = -1;
% Construct initial geometry
for j=1:nely
    for i=1:nelx-mod(j+1,2)
        xi((j-1)*nelx+i-floor((j-1)/2)) = Diam*(i-0.5)+(1-mod(j,2))*Diam/2;
        yi((j-1)*nelx+i-floor((j-1)/2)) = Diam*((j-1)*cos(pi/6)+0.0);
    end
end
% compute neighbour list and equilibrium distances
[N,nn,Li]=nlist(xi,yi,cut);
% initialize variables and plotting
loop = 0; 
change = 1.;
x = xi; y = yi;
% cc=[1-m,1-m,1-m];
% scatter(x, y, Diam*100, cc, 'filled')
% linkdata on
% START ITERATION
 while change > 0.001
    loop=loop+1;
    mold = m;
    % DEM-ANALYSIS assuming harmonic potential
    [x,y]=DEM(m,x,y,nn,N,nelx,nely,kspr,Li,Fxe,Fye,damp,Diam);
    % compute strain energy (cost function) and sensitivity
    Eer=0;
    dc = zeros(length(xi),1);
    for i=1:length(xi)
        for s=1:nn(i)
            L = sqrt( (x(i)-x(N(i,s)))^2 + (y(i)-y(N(i,s)))^2 );
            Eer = Eer + 1/4 * m(i)^2*m(N(i,s))^2 * kspr * ( L - Li(i,s) )^2; % 1/4 cause double counting each interaction
            %dc(i) = dc(i) - 1/4 * m(N(i,j))*(2*m(i)+m(N(i,j)))  * kspr * ( L - Li(i,j) )^2;
            %dc(i) = dc(i) - 1/2 * m(N(i,j))  * kspr * ( L - Li(i,j) )^2;
            dc(i) = dc(i) - m(N(i,s))^2*m(i)  * kspr * ( L - Li(i,s) )^2;
        end
    end
%     fileID = fopen('dump.txt','w')
    % NOTE: filtering of sensitivities (aka coarse graining) should go here...
    %[dc]   = check(nelx,nely,rmin,rho,dc); 
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [m]    = OC(m,mass,dc);
    %PRINT RESULTS
    change = max(abs(m-mold));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',Eer) ...
       ' Mass.: ' sprintf('%6.3f',sum(m)/(length(xi))) ...
        ' ch.: ' sprintf('%6.3f',change )])
    % PLOT DENSITIES  
    %colormap(gray); imagesc(-m); axis equal; axis tight; axis off;pause(1e-6);
       cc=[1-m,1-m,1-m];
       linkdata on
       scatter(x, y, Diam*100, cc, 'filled')  
 end
 % Write Dumpfiles
 fileID = fopen('outputs.txt','w');
 fprintf(fileID,'%6s %12s\n','Maxdisp','Eer');
 fprintf(fileID,'%6.6f %12.6f\n',abs(y(ceil(nelx/2))),Eer);
 fclose(fileID);
 toc

%  Eer=zeros(length(nonzeros(N)/2),1);
%  d=0;
%      for i=1:length(xi)
%         for s=1:nn(i)
%             j= N(i,s);
%             if j>i && j~=0
%                 d=d+1;
%                 L = sqrt( (x(i)-x(j))^2 + (y(i)-y(j))^2 );
%                 Eer(d) = 1/2 * m(i)^2*m(N(j))^2 * kspr * ( L - Li(i,s) )^2;
%             end
%         end
%      end
%      Eer = Eer./max(Eer)
%      d=0;
%      for i=1:length(xi)
%         for s=1:nn(i)
%             j= N(i,s);
%             if j>i && j~=0
%                 d=d+1;
%                 linkdata on
%                 plot([x(i);x(j)],[y(i);y(j)],'color',[Eer(d),1-Eer(d),0],'linewidth',4*m(i)*m(j)+0.001)
%                 hold on
%             end
%         end
%      end
     

% Building neighbour list and computing initial equilibrium lengths 
function [N,nn,Li]=nlist(x,y,cut)
nn=zeros(length(x),1);
for i=1:length(x)
    for j=1:length(x)
        L = sqrt( (y(i)-y(j))^2 + (x(i)-x(j))^2 );
        if i~=j && L < cut
            nn(i)=nn(i)+1;
            N(i,nn(i))=j;
            Li(i,nn(i))=L;
        end
    end
end
% DEM energy minimisation by steepest descent
function [x,y]=DEM(m,x,y,nn,N,nelx,nely,kspr,Li,Fxe,Fye,damp,Diam);
tol = Diam;
count = 0;
for i=1:length(m)
    for s=1:nn(i)
        j = N(i,s);
        k(i,s) = m(i)^2 *m(j)^2  * kspr;
    end
end
while tol > Diam*nely/1e8
    count = count + 1;
    tol = 0; 
    Fx=Fxe;  Fy = Fye;
    for i=1:length(m)
        for s=1:nn(i)
            j = N(i,s);
            if j > i
                dx = x(j)-x(i);     dy = y(j)-y(i);
                L = sqrt( dx^2 + dy^2 ); 
                F = k(i,s) * (L - Li(i,s));
                Fxx = F * dx/L;     Fyy = F * dy/L;
                Fx(i) = Fx(i) + Fxx;    Fx(j) = Fx(j) - Fxx;
                Fy(i) = Fy(i) + Fyy;    Fy(j) = Fy(j) - Fyy;
            end
        end
    end
    %constraints and displacement
    Fx(1)=0;
    Fy(1)=0;
    Fy(nelx) = 0;
%     for i=1:nely
%         pos = 1+nelx*(i-1);
%         Fx(pos) = 0;
%     end
    x = x + Fx/damp;
    y = y + Fy/damp;
    %constraints and criterion
    tol = max(abs(Fx)+abs(Fy))/damp;  
    % interim output to appreciate speed of code
     if mod(count,10000)==0
%         disp(sprintf('2 %f %e %e',y(1),Fy(1),tol));
    end
end
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mnew]=OC(m,mass,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  mnew = max(0,max(m-move,min(1.,min(m+move,m.*sqrt(-dc./lmid)))));
  if sum(mnew) - mass * length(m) > 0
    l1 = lmid;
  else
    l2 = lmid;
  end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
  for j = 1:nely
    sum=0.0; 
    for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
      for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
        fac = rmin-sqrt((i-k)^2+(j-l)^2);
        sum = sum+max(0,fac);
        dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
      end
    end
    dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
  end
end


    
    
%
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
