close all
clear all
clc

% dt = 2.5e-11;
% endt = 12e-9;
% t = 0:dt:endt;
% n = length(t);
% 
% q = 0.05;
% fc = 1*1e9;
% 
% fs = 1/dt / 1e3 * 2;
% f = fs *(0:n/2 + 1);
% 
% v = rickerPulse(t,q,fc);
% 
% Vf = fft(v);
% magVf = abs(Vf);
% p1 = magVf(1:n/2+2);
% p1(2:end-1) = 2*p1(2:end-1);
% 
% figure()
% plot(t,v);
% 
% figure()
% plot(f, p1)

figure()
i = 1;
f1 = 1;
f2 = 2;
fn = f2 * 2;
Tn = 1/fn;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
nsteps = 41;
tend = 5;

t = linspace(0,tend,nsteps);

x = sin(w1*t) + sin(w2*t);
Xf = fft(x,nsteps);
magXf = abs(Xf);

Fs = nsteps / tend;
n = length(Xf);
len = n+1;

sampleRate = n / t(n);
Ts = 1/sampleRate;

fs = sampleRate * (-(len-1)/2:(len-2)/2)/len;
Yf = zeros(1,n);
Yf(1:n/2) = Xf(n/2+1:end);
Yf(n/2+1:end) = Xf(1:n/2);
magYf = abs(Yf);

%Yf = reArrange(Xf);

subplot(3,3,i)
i = i+1;
plot(t,x)

subplot(3,3,i)
i = i+1;
plot(abs(Xf))

subplot(3,3,i)
i = i+1;
plot(fs,magYf)