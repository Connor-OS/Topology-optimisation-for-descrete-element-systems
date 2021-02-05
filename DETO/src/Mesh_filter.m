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
dc = dcn;