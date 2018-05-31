function [alpha, beta] = calcGamma(eps, f)
omega = 2 * pi * f;
u0 = 4*pi*1e-7;
epsRe = real(eps);
epsIm = imag(eps);
lossTangent = epsIm ./ epsRe;

alpha = omega * sqrt(u0*epsRe/2.*(sqrt(1+lossTangent.^2)-1));
beta = omega * sqrt(u0*epsRe/2.*(sqrt(1+lossTangent.^2)+1));
