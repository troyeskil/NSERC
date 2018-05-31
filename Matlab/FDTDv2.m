close all
clear all
clc

%Constants
eps0 = 8.854e-12;
u0 = 4*pi*1e-7;
c0 = 1/sqrt(eps0*u0);
eta0 = sqrt(u0 / eps0);

%Wave Properties
f = 1e9; %[Hz]
Emag = 1; %[V/m]
EmagNorm = Emag / eta0;

omega = 2*pi*f;
k = omega / c0;
lambda0 = 2 * pi / k;

%Simulation Parameters
d = 0.5; %height of simulation
nsteps = 400;

dz = lambda0 / 20/2;
dt = dz / c0 / 2;
zs = 0:dz:d;
ncells = length(zs);

S = c0 * dt / dz; %Courant number
Nlambda = lambda0 / dz; % should be greater than Ntransition
Ntransition = 2 * pi * S / (acos(1 - 2 * S.^2));
theta = 1 + (1/S)^2 * (cos(2*pi*S / Nlambda) - 1);

%Dielectric
epsr = 4;
d1 = d/2;
dielCond = 0.05;

%Numerical Wave
k_num = 1 / dz * acos(theta);
vp_num = omega ./ real(k_num);

%Material Properties Mesh
eps = ones(1,ncells) * eps0;
eps(ceil(ncells/2):ncells) = eps0*epsr;
conds = zeros(1,ncells);
conds(ceil(ncells/2):ncells) = dielCond;


%Initialize Fields
Ex = zeros(1, ncells);
Hy = zeros(1, ncells);


M = moviein(nsteps);
t = 0:dt:dt*(nsteps - 1);

cc = c0*dt/dz;

%Set up ABC
oldE2 = Ex(ncells-1);
oldE1 = oldE2;

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
    Ex(1) = EmagNorm * sin(omega * t(i));
    Ex(2:ncells-1) = k1(2:ncells-1) .* Ex(2:ncells-1) + k2(2:ncells-1) .* (Hy(1:ncells-2) - Hy(2:ncells-1));
    
    %Absorbing Boundary Conditions
    Ex(ncells) = oldE2;
    oldE2 = oldE1;
    oldE1 = Ex(ncells - 1);
    
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

   
