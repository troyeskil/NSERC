function [Yf,fs] = reArrange(Xf,t)

n = length(t);
len = n+1;

sampleRate = n / t(n);
Ts = 1/sampleRate;
fs = sampleRate * (-(len-1)/2:(len-2)/2)/len;

Yf = zeros(1,n);
Yf(1:n/2) = Xf(n/2+1:end);
Yf(n/2+1:end) = Xf(1:n/2);
end