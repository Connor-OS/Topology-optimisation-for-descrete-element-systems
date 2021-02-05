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