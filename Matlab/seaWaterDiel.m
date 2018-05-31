function [epsr] = seaWaterDiel(f,Tk,S)
%Pg 125

eps0 = 8.854e-12;
Tc = Tk - 273.15;
%S = S/100; %PSU
Spsu = S * 100;
%Tc = 0;

%Tabel of Coefficients
a1 = 0.46606917e-2;
a2 = -0.26087876e-4;
a3 = -0.63926782e-5;
a4 = 0.63000075e1;
a5 = 0.26242021e-2;
a6 = -0.42984155e-2;
a7 = 0.34414691e-4;
a8 = 0.17667420e-3;
a9 = -0.20491560e-6;
a10 = 0.58366888e3;
a11 = 0.12684992e3;
a12 = 0.12684992e3;
a13 = 0.38957681e-6;
a14 = 0.30742330e3;
a15 = 0.12634992e3;
a16 = 0.37245044e1;
a17 = 0.92609781e-2;
a18 = -0.26093754e-1;

%Relaxation frequencies @273k
% f01 = 8.9e9;
% f02 = 201.8e9;

epsw0 = 87.85306*exp(-0.00456992*Tc - a1*S - a2*S^2 - a3*S*Tc);
epsw1 = a4*exp(-a5*Tc - a6*S - a7*S*Tc);

tauw1 = (a8 + a9*S)*exp(a10/(Tc + a11))*1e-9;
tauw2 = (a12+a13*S)*exp(a14/(Tc + a15))*1e-9;

epswInf = a16 + a17*Tc + a18*S;

%Conductivity
sigmaT = 2.903602 + 8.607*10^-2*Tc + 4.738817*10^-4*Tc^2 - 2.991*10^-6*Tc^3 + 4.3041*10^-9*Tc^4;
pS = Spsu*(37.5109 + 5.45216*Spsu - 0.014409*Spsu^2)/(1004.75 + 182.283*Spsu + Spsu^2);

alpha0 = (6.9431 + 3.2841*Spsu - 0.099486*Spsu^2)/(84.85 + 69.024*Spsu + Spsu^2);
alpha1 = 49.843 - 0.2276*Spsu + 0.00198*Spsu^2;
qTS = 1 + alpha0 * (Tc - 15)/(Tc + alpha1);

sigmai = sigmaT*pS*qTS;

epsrRe = epswInf + (epsw0-epsw1)/(1+(2*pi*f*tauw1)^2) + (epsw1-epswInf)/(1+(2*pi*f*tauw2)^2);
epsrIm = 2*pi*f*tauw1*(epsw0-epsw1)/(1+(2*pi*f*tauw1)^2) + 2*pi*f*tauw2*(epsw1-epswInf)/(1+(2*pi*f*tauw2)^2) + sigmai/(2*pi*eps0*f);

epsr = complex(epsrRe,epsrIm);
end