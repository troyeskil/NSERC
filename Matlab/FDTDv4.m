close all
clear all
clc

%Constants
eps0 = 8.854e-12;
u0 = 4*pi*1e-7;
c0 = 1/sqrt(eps0*u0);
eta0 = sqrt(u0 / eps0);

%Simulation Inputs
speedFactor = 1; %Only display every nth frame
fps = 48; %default 12
d = 0.5; %[m] height of simulation
f = 1e9; %[Hz]
Tk = 273.15-10; %Temperature, Kelvin
snowDensity = 0.3;%[g/cm^3] Typical value 0.1 < x < 0.4.  Taken from figure 4-17, pg 141.

depths = [0.2, 0.3];
diels = ["snow", "ice"];

%Wave properties
omega = 2*pi*f;
Emag = 1; %[V/m]
EmagNorm = Emag / eta0;
k = omega / c0;
lambda0 = 2 * pi / k;

nlayers = length(diels);
[epsLayers,~] = generateGrids(nlayers, nlayers + 1, [1:nlayers], diels, Tk, f, snowDensity);
epsrLayers = epsLayers / eps0;
epsrMax = max(epsrLayers);
lambdaMin = lambda0 / sqrt(epsrMax);

%Simulation Parameters
dz = lambdaMin / 20;
dt = dz / c0 / 2;
zs = 0:dz:d;
ncells = length(zs);
nsteps = ncells * 5;
t = 0:dt:dt*(nsteps - 1);

[eps, conds] = generateGrids(d, ncells, depths, diels, Tk, f, snowDensity);

%Source
q = 0.05; %damping factor, 0.2 for 250MHz, 0.05 for 1GHz
Esource = Emag * sinPulse(t,f);
%Esource = Emag * rickerPulse(t,q,f);
%Esource = Emag * sin(omega*t);

%Initialize Fields
Ex = zeros(1, ncells);
Hy = zeros(1, ncells);


M = moviein(floor(nsteps/speedFactor));

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

%Reciever
echo = zeros(1,nsteps);

%Simulation
p = simPlots(Ex, depths, diels, zs, d, ncells, nlayers);
j = 0;

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
    
    %Recievd Field
    echo(i) = Ex(1) * eta0;
    
    j = j + 1;
    
    %Plots
    if (j >= speedFactor)
        j = 0;
        ex = Ex * eta0;
        p.YData = ex;
        M(:,round(i/speedFactor)) = getframe;
    end

end

movie(M,1,fps);

csvwrite('echoSignal.csv',echo);
csvwrite('timeStamps.csv',t);

%Plot recieved signal
figure()
plot(t,echo)

function [eps, condEff] = generateGrids(d, ncells, depths, diels, Tk, f, snowDensity)
eps0 = 8.854e-12;
epsr = ones(1,ncells);
condEff = zeros(1,ncells);
for i = 1:length(depths)
    layer = diels(i);
    if (layer == "snow")
        epsLayer = drySnowDiel(snowDensity, Tk, f);
    elseif (layer == "ice")
        epsLayer = pureIceDiel(Tk, f);
    end
    startCell = mapToGrid(d,ncells, depths(i));
    epsr(startCell:ncells) = epsLayer;
end
eps = real(epsr) * eps0;
condEff = 2 * pi * f * imag(epsr) * eps0;
end

function [vertBounds, diel] = getPlotParams(d, ncells, depths, diels, zs)
h = 2;
vertBounds = [-h, h];
loc = mapToGrid(d, ncells, depths);

labels = string(zeros(1,nlayers+1));
labels(nlayers+1) = "Efield";
dielXs = zeros(1,nlayers*2);
dielYs = zeros(1,nlayers*2);

for i = 1:nlayers
    z = [zs(loc(i)),zs(loc(i))];
    dielXs(2*i-1:2*i) = z;
    dielYs(2*i-1:2*i) = vertBounds;
    labels(i) = diels(i);
end
plot(dielXs,dielYs,'--')
end

function [p] = simPlots(ex, depths, diels, zs, d, ncells, nlayers)

h = 2;
vertBounds = [-h, h];
loc = mapToGrid(d, ncells, depths);

labels = string(zeros(1,nlayers+1));
labels(1) = "Efield";

figure(1)
p = plot(zs, ex);
hold on

for i = 1:nlayers
    z = [zs(loc(i)),zs(loc(i))];
    plot(z, vertBounds,'--')
    hold on
    labels(i+1) = diels(i);
end

legend(labels)
hold on
end





function [signal] = sinPulse(t, f)
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