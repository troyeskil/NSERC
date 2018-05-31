function [epsr] = brineDiel(f, Tk)

eps0 = 8.854e-12;
Tc = Tk - 273.15;
delta = 25 - Tc;
Si = 0.04; %in ppt

if Tc <= -2
    if Tc >= -8.2
        Sb = 1.725 - 18.756*Tc - 0.3964*Tc^2;
    elseif Tc >= -22.9
        Sb = 57.041 - 9.929*Tc - 0.16204*Tc^2 - 0.002396*Tc^3;
    elseif TC >= -36.8
        Sb = 242.92 + 1.5299*Tc + 0.0429*Tc^2;
    elseif Tc >= -43.2
        Sb = 508.18 + 14.535*Tc + 0.2018*Tc^2;
    end
end

vb = 1e-3*Si*(-49.185/Tc + 0.532);

epswInf = 4.9;

Nb = Sb*(1.707e-2 + 1.205e-5*Sb + 4.058e-9*Sb^2);

epsBT0 = 88.045-0.4147*Tc + 6.295*1e-4*Tc^2 + 1.075*1e-5*Tc^2;
a1 = 1 - 0.255*Nb + 5.15e-2*Nb^2 - 6.89e-3*Nb^3;

ktaub0 = 1.1109*1e-10 - 3.824*1e-12*Tc + 6.938*1e-14*Tc^2 - 5.096*1e-16*Tc^3;
taub0 = ktaub0 / 2 / pi;
b1 = 1 + 0.146e-2*Nb*Tc - 4.89e-2*Nb - 2.97e-2*Nb^2 + 5.64e-3*Nb^3;

sigmab25 = Nb * (10.39 - 2.378*Nb + 0.683*Nb^2 - 0.135*Nb^3 + 1.01e-2*Nb^4);
c1 = 1 - 1.96e-2*delta + 8.08e-5*delta^2 - Nb*delta*(3.02e-5 + 3.92e-5*delta + Nb*(1.72e-5 - 6.58e-6*delta));

epsb0 = epsBT0 * a1;
taub = taub0 * b1;
sigmab = sigmab25 * c1;

epsrRe = epswInf + (epsb0 - epswInf)/(1 + (2*pi*f*taub)^2);
epsrIm = 2*pi*f*taub*(epsb0-epswInf)/(1+(2*pi*f*taub)^2) + sigmab/(2*pi*eps0*f);

epsr = complex(epsrRe,epsrIm);
        