damp = max(300,3*kspr);
tol = Diam;
count = 0;
L=Li;
while tol > Diam*nely/(kspr*2e7)
    count = count + 1;
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