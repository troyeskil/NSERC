function [eps] = waterDiel(f, Tk)

eps0 = 8.854e-12;
Tc = Tk - 273.15;

epsw0 = 88.045-0.4147*Tc + 6.295*1e-4*Tc^2 + 1.075*1e-5*Tc^2;
epswInf = 4.9;

ktau = 1.1109*1e-10 - 3.824*1e-12*Tc + 6.938*1e-14*Tc^2 - 5.096*1e-16*Tc^3;

epsrRe = epswInf + (epsw0-epswInf)/(1+(ktau)^2);
epsrIm = f*ktau*(epsw0-epswInf)/(1+(f*ktau)^2);

eps = complex(epsrRe,epsrIm) * eps0;