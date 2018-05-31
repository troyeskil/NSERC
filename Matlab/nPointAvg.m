function xbar = nPointAvg(x, n)

k = ceil(n/2);
len = length(x);
xbar = zeros(1,len);

for j = 1:n
    xbar(k:len-k+1) = xbar(k:len-k+1) + x(j:len-n+j);
    xbar(j) = sum(x(1:j))/j;
    xbar(len-j+1) = sum(x(len-j+1))/j;
end
