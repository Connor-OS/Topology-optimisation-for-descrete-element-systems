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
m = mnew;