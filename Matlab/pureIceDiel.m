function [eps] = pureIceDiel(Tk,f)
epsReal = 3.1884; %+9.1e-4*T

freq = f / 1e9;
theta = 300/Tk - 1; %T in K
B1 = 0.0207; %K/GHz
b = 335; %[K]
B2 = 1.16e-11; %[GHz^-3]
alpha0 = (0.00504+0.0062*theta) * exp(-22.1*theta);
beta0 = B1/Tk*exp(b/Tk)/(exp(b/Tk)-1)^2 + B2 * freq^2 + exp(-9.963+0.0372*(Tk - 273.16));
epsIm = alpha0 / freq + beta0 * freq;

eps = epsReal + 1j*epsIm;
end