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
Qmin
stress_perpart
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