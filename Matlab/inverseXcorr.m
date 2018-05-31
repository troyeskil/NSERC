function [a] = inverseXcorr(b,r)
n = length(b);
corrLength = length(r);
Bf = fft(b,corrLength);
Rf = fft(r,corrLength);
Af = Rf .* conj(Bf);
a = (ifft(Af));
%a = a(n:end);

% a = deconv(b,r);