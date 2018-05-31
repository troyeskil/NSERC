% This is a 1D FDTD simulation with pulse
% It displays a "movie" of the signal
% Size of the FDTD space
clear;
ke=50;
% Position of the source
ks=1;
% Number of time steps
nsteps=100;
% Cell size and time stepping
c0=3.e8;
dx=0.01;
dt=dx/(2.*c0);
% Constants
cc=c0*dt/dx;
% Initialize vectors
ex=zeros(1,ke);
hy=zeros(1,ke);
% Gaussian pulse
t0=20;
spread=8;
% Start loop
M=moviein(nsteps);
for t=1:nsteps
    % E field loop
    for k=2:ke-1
        ex(k)=ex(k)+cc*(hy(k-1)-hy(k));
    end
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