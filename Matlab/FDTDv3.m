close all
clear all
clc

%Constants
eps0 = 8.854e-12;
u0 = 4*pi*1e-7;
c0 = 1/sqrt(eps0*u0);
eta0 = sqrt(u0 / eps0);

%Simulation Inputs
d = 0.5; %height of simulation
nsteps = 400;
f = 1e9; %[Hz]
Tk = 273.15-10; %Temperature, Kelvin

%Wave properties
omega = 2*pi*f;
Emag = 1; %[V/m]
EmagNorm = Emag / eta0;
k = omega / c0;
lambda0 = 2 * pi / k;

%Simulation Parameters
dz = lambda0 / 20/2;
dt = dz / c0 / 2;
zs = 0:dz:d;
ncells = length(zs);
t = 0:dt:dt*(nsteps - 1);

%Dielectric
epsr = 16;
d1 = d/2;
dielCond = 0;

%Source
q = 0.05; %damping factor, 0.2 for 250MHz, 0.05 for 1GHz
Esource = Emag * sinPulse(t,f);
%Esource = Emag * rickerPulse(t,q,f);
%Esource = Emag * sin(omega*t);

%Snow
snowDensity = 0.3;%[g/cm^3] Typical value 0.1 < x < 0.4.  Taken from figure 4-17, pg 141.
epsIce = pureIceDiel(Tk,f);
epsDrySnow = drySnowDiel(snowDensity, epsIce);

%Material Properties Grid
eps = ones(1,ncells) * eps0;
eps(ceil(ncells/2):ncells) = eps0*epsr;
conds = zeros(1,ncells);
conds(ceil(ncells/2):ncells) = dielCond;


%Initialize Fields
Ex = zeros(1, ncells);
Hy = zeros(1, ncells);


M = moviein(nsteps);

%Unused code
%cc = c0*dt/dz;
%
% k_num = 1 / dz * acos(theta);
% vp_num = omega ./ real(k_num);
% 
% S = c0 * dt / dz; %Courant number
% Nlambda = lambda0 / dz; % should be greater than Ntransition
% Ntransition = 2 * pi * S / (acos(1 - 2 * S.^2));
% theta = 1 + (1/S)^2 * (cos(2*pi*S / Nlambda) - 1);

%Set up ABC
%Far Side
oldE2far = Ex(ncells-1);
oldE1far = oldE2far;
%Close Side
oldE2LHS = Ex(1);
oldE1LHS = oldE2LHS;
oldSource1 = 0;
oldSource2 = oldSource1;

%Constant Terms used in calculations
epsr = eps/eps0;
q1 = conds * dt ./ (2 * eps);
q2 = c0 * dt / dz;

k1 = (1 - q1) ./ (1 + q1);
k2 = q2 ./ epsr ./ (1 + q1);
k3 = 1;
k4 = q2;

%Simulation
for i = 1:(nsteps)
    
    %Efield calculations
    Ex(1) = EmagNorm * Esource(i) + oldE2LHS;
    Ex(2:ncells-1) = k1(2:ncells-1) .* Ex(2:ncells-1) + k2(2:ncells-1) .* (Hy(1:ncells-2) - Hy(2:ncells-1));
    
    %Absorbing Boundary Conditions
    %Far Side
    Ex(ncells) = oldE2far;
    oldE2far = oldE1far;
    oldE1far = Ex(ncells - 1);
    
    %Close Side
    oldE2LHS = oldE1LHS;
    oldE1LHS = Ex(2) - oldSource2;
    oldSource2 = oldSource1;
    oldSource1 = EmagNorm * Esource(i);
    
    %H Field calculations
    Hy(1:ncells-1) = k3 .* Hy(1:ncells-1) + k4 .* (Ex(1:ncells-1) - Ex(2:ncells));
    
    ex = Ex * eta0;
    plot(ex); axis([1 ncells -2 2]);
    M(:,i) = getframe;
    %input('')

end

movie(M,1);


function [i] = mapToGrid(gridLength, gridPoints, loc)

pointsPerMeter = gridPoints / gridLength;
i = ceil(pointsPerMeter * loc);
end

function [eps] = pureIceDiel(Temp,f)
epsReal = 3.1884; %+9.1e-4*T

freq = f / 1e9;
theta = 300/Temp - 1; %T in K
B1 = 0.0207; %K/GHz
b = 335; %[K]
B2 = 1.16e-11; %[GHz^-3]
alpha0 = (0.00504+0.0062*theta) * exp(-22.1*theta);
beta0 = B1/Temp*exp(b/Temp)/(exp(b/Temp)-1)^2 + B2 * freq^2 + exp(-9.963+0.0372*(Temp - 273.16));
epsIm = alpha0 / freq + beta0 * freq;

eps = epsReal + j*epsIm;
end

function [eps] = drySnowDiel(snowDensity, epsIce)

%epsrIce = 3.17;%Ignoring the imaginary part
vi = snowDensity / 0.9167;
eps = 1 + 3*vi*(epsIce - 1) / ((2+epsIce) - vi*(epsIce-1));
end

function [signal] = sinPulse(t, f);
    period = 1/f;
    maskSinPulse = t < period;
    signal = sin(2*pi*f*t) .* maskSinPulse;
end

function [signal] = rickerPulse(t, q, fc)

T = (0.8 + (1-q)/15)/fc;

vt0 = vPulse(t, T);
vt1 = vPulse(t - T/2, T);
vt2 = vPulse(t - T, T);

signal = vt0 - (2-q) * vt1 + (1-q) * vt2;
end

function [pulse] = vPulse(t, T)

mask = t < T;
pulse = 1/2*(1+cos(pi*(t-T/2)/(T/2)));
pulse = pulse .* mask;
end