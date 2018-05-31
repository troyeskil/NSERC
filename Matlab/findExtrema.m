function loc = findExtrema(lo,hi,x)
xprime = diff(x);
%xprime = nPointAvg(xprime,7);

if lo >= hi - 1
    loc = lo;
else
    if xprime(lo) > 0
        mask = xprime < 0;
        loc = lo - 1 + find(mask(lo:hi-1),1);
    else
        mask = xprime > 0;
        loc = lo - 1 + find(mask(lo:hi-1),1);
    end
end