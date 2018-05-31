function [f] = inverseConv(g,y)
n = length(y);
Yf = fft(y);
Gf = fft(g,n);
Ff = Yf ./ Gf;
f = ifft(Ff);