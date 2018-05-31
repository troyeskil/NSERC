% This is a 1D FDTD simulation with pulse
% It displays a "movie" of the signal
% Size of the FDTD space
close all
clear all
clc

ke=50;
% Position of the source
ks=1;
% Number of time steps
nsteps=400;

%Constants
eps0 = 8.854e-12;
u0 = 4*pi*1e-7;
c0=3.e8;
eta0 = sqrt(u0/eps0);

%Dielectric Boundary
kloc = 25;
epsDiel = 16;
dielcond = 0;
vp = zeros(1,ke) + c0;
vp(kloc:ke) = c0 / sqrt(epsDiel);

%Material Properties
conds = zeros(1,ke);
conds(kloc:ke) = dielcond;

eps = ones(1,ke) * eps0;
eps(kloc:ke) = eps0 * epsDiel;
epsr = eps / eps0;

% Cell size and time stepping
dx=0.005;
dt=dx/(2.*c0);

% Initialize vectors
ex=zeros(1,ke);
hy=zeros(1,ke);

%Sin wave
ts = 0:dt:dt*nsteps;
f = 1e6;
eSin = sin(2*pi*f*ts)/eta0;

% Gaussian pulse
t0=20;
spread=8;
ePulse = exp(-.5*((ts)/spread).^2);

%ABS Conditions
exn1 = ex(ke);
exn2 = exn1;

% Constants
q1 = conds * dt ./ (2 * eps);

cc = c0*dt/dx;
ca = (1 - q1)./(1 + q1);
cb = (cc ./ epsr) ./ (1 + q1);

% Start loop
M=moviein(nsteps);
for t=1:nsteps
    % E field loop
    for k=2:ke-1
        ex(k)= ca(k)*ex(k)+cb(k)*(hy(k-1)-hy(k));
    end
    %Absorbing Boundary at right side
    ex(ke) = exn2;
    exn2 = exn1;
    exn1 = ex(ke-1);
    
    % Source
    ex(ks)=exp(-.5*((t-t0)/spread)^2);
    % H field loop
    for k=1:ke-1
        hy(k)=hy(k)+cc*(ex(k)-ex(k+1));
    end
    plot(ex);axis([1 ke -2 2]);
    M(:,t) = getframe  ;
    %  input('')
end
movie(M,1);