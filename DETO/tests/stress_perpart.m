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
        W(i,6) = sqrt(W(i,1)^2+W(i,2)^2+3*W(i,3)^2-W(i,1)*W(i,2)); % Deviatoric stress
    end
end
pvol =((pi*Diam^2)/4)/0.9069;
W = (W/pvol);