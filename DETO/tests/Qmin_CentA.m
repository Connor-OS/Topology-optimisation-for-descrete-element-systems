vx = zeros(length(x),1); vy = vx;
freeze_step = 0;
tol_min = 5e-10;
tol = 2*tol_min;
nstep = 0;
dmax = 0.01;
dt = 0.01;
Eer = 0;
low_m = min(m);
% Euler integration
while tol > tol_min
    nstep = nstep + 1;
    Fx = Fxe;  Fy = Fye;
    % Compute force
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
    % Constraints
    Fy(ceil(length(m)/2)-floor(nelx/2)) = 0;
    Fy(ceil(length(m)/2)+floor(nelx/2)) = 0;
    % Compute scale factor
    vdotf = dot(vx,Fx) + dot(vy,Fy);
    % Overshoot check
    if vdotf < 0
        freeze_step = nstep;
        vx = zeros(length(x),1); vy = vx;
    else
        % Scale velocities
        fdotf = dot(Fx,Fx) + dot(Fy,Fy);
        if fdotf == 0; scale = 0;
        else scale = vdotf/fdotf;
        end
        vx = scale .* Fx;
        vy = scale .* Fy;
    end
    % Limit dmax
    vmax = max(max(abs(vx)),max(abs(vy)));
    if dt*vmax > dmax; dt = dmax/vmax;
    end
    % Euler integration step
    x = x + dt .* vx;
    y = y + dt .* vy;
    for i = 1:length(x)
        if m(i) > 0
            vx(i)= vx(i)+dt/m(i)*Fx(i);
            vy(i)= vy(i)+dt/m(i)*Fy(i);
        else
            vx(i) = 0;
            vy(i) = 0;
        end
    end
    % Tolerance
    Eold = Eer;
    Eer = 0;
    for i = 1:length(xi)
        for s = 1:nn(i)
            j = N(i,s);
            Eer = Eer+1/4*m(i)^2*m(j)^2*kspr*(L(i,s)-Li(i,s))^2;
        end
    end
    if nstep - freeze_step < 5
        tol = 2*tol_min;
    else
        tol = abs(Eer-Eold)/Eold;
    end
end