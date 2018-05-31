close all
clear all
clc

f = 20e9;
Tk = 273.15 + 20;
Tb = 273.15 - 5;
S = 0.3254;
eps0 = 8.854e-12;

epsr1 = seaWaterDiel(f,Tk,S) / eps0;
epsr2 = seaWaterDiel(f,Tk,0) / eps0;
epsr3 = waterDiel(f, Tk) / eps0;

n = 100;
flow = 9;
fhi = 12;
fs = logspace(flow,fhi,n);
fs2 = logspace(9,11,n);

epsS = ones(1,n) * 1i;
epsW = ones(1,n) * 1i;
epsB = ones(1,n) * 1i;

for i = 1:n
    freq = fs(i);
    epsS(i) = seaWaterDiel(freq, Tk, S) / eps0;
    epsW(i) = seaWaterDiel(freq, Tk, 0) / eps0;
    epsB(i) = brineDiel(fs2(i),Tb) / eps0;
end

figure()
subplot(2,2,1)
loglog(fs, real(epsS))
hold on
loglog(fs, real(epsW))

subplot(2,2,2)
loglog(fs, imag(epsS))
hold on
loglog(fs, imag(epsW))

subplot(2,2,3)
loglog(fs2, real(epsB))
hold on
loglog(fs2, imag(epsB))