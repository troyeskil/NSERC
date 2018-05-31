%function [diels, depths] = processEcho(t, rawSignal, f, sentPulse)

close all
clear all
clc

i = 1;
f = 10e9;
d = 2.5;
ncells = 100;

eps0 = 8.854e-12;
u0 = 4*pi*1e-7;
c0 = 1/sqrt(eps0*u0);

%imports
rawSignal = csvread('echoSignal.csv');
t = csvread('timeStamps.csv');
sentPulse = csvread('sentPulse.csv');

figure(1)
subplot(2,2,i)
i = i+1;
plot(t,rawSignal);
title("Signal Recieved")

subplot(2,2,i)
i = i+1;
plot(t,rawSignal - sentPulse);
title("Echo Recieved")

%Initial Calculations
pulseT = 1/f;
nsteps = length(t);

%Processing of Signal
locPulse = findExtrema(1,nsteps,sentPulse);
shiftPulse = locPulse - 1;

signal = rawSignal - sentPulse;
keepLoop = 1;
k = 2;

minPk = 5;
shift = 150;
[newRmm, lags] = xcorr(signal, sentPulse);
[rPulse, ~] = xcorr(sentPulse, sentPulse);

while keepLoop
    
    [absPks,locs] = findpeaks(abs(newRmm),'MinPeakHeight', minPk);
    
    if (k > 5 || length(locs) < 5)
        keepLoop = 0;
    else
        
        %%%%%
        %Echo removal using cross correlation
        %%%%%
        
        %Peaks of signal
        Rmm = newRmm;
        lengthcorr = length(Rmm);
        pks = Rmm(locs);
        
        %Peaks of sent pulse
        [pulseAbsPks, pulseLocs] = findpeaks(abs(rPulse), 'MinPeakHeight', minPk);
        pulsePks = rPulse(pulseLocs);
        
        mid = find(lags>0,1);
        scale = pks(1) / pulsePks(1);
        
        scales(k-1) = scale;
        echoLocations(k-1) = locs(1);
        
        rEcho = zeros(1,lengthcorr);
        rEcho(locs(1)-shift:end) = rPulse(pulseLocs(1)-shift:-locs(1)+lengthcorr+pulseLocs(1)) * scale;
        
        newRmm = Rmm - rEcho;

        %keepLoop = 0;
        k = k+1;
        j = 1;
        
        %%%%%
        %Plotting
        %%%%%
        figure()
        
        subplot(2,2,j)
        j = j+1;
        plot(Rmm)
        title("Rmm")
        hold on
        scatter(locs,pks,'ko')
        hold off
        
        subplot(2,2,j)
        j = j+1;
        plot(rPulse)
        title("rPulse")
        hold on
        scatter(pulseLocs,pulsePks,'ko')
        hold off
        
        subplot(2,2,j)
        j = j+1;
        plot(rEcho)
        title("rEcho")
        
        subplot(2,2,j)
        j = j+1;
        plot(newRmm)
        title("newRmm")
        
    end
end

sentPower = sum(sentPulse.^2);
echoPower = sum(signal.^2);
nrcs = echoPower / sentPower;

% figure(1)
% thresholdLine = ones(1,corrLength)*thresh;
% subplot(2,2,i)
% i = i+1;
% plot(r)
% hold on
% plot(thresholdLine)
% hold on
% plot(-thresholdLine)
% hold off


start = find(lags>0);
locations = echoLocations - pulseLocs(1) + 1;
times = t(locations);
nlayers = length(scales);
errorRelative = zeros(1,nlayers);
reflectCoefs = scales;

alpha = zeros(1,nlayers);
beta = zeros(1,nlayers);

layers = string(zeros(1,nlayers));
dielEpsr = ones(1,nlayers+1);

Tk = 273.15 - 10;
snowDensity = 0.3;

dielTypes = ["snow","ice","air","sea water","water", "brine","lossless 4"];
Tks = ones(1,length(dielTypes)) * Tk;
nTypes = length(dielTypes);
epsrSpecified = ones(1,nTypes);
epsrSpecified(end) = 4;
[epsRe, condEff, epsIm, ~] = generateGrids(nTypes + 1, nTypes + 1, (1:nTypes) - 0.5, dielTypes, epsrSpecified, Tks, f, snowDensity);
epsrTypes = complex(epsRe,epsIm) / eps0;
epsrTypes = epsrTypes(2:end);

for i = 1:length(scales)
    
    if i>1
        %Account for 2 way transit of wave
        gamma2 = calcReflect(dielEpsr(i),dielEpsr(i-1));
        gamma2back = abs(calcReflect(dielEpsr(i-1),dielEpsr(i)));
        atten = real(exp(-alpha(i-1)*2*depths(i-1)));
        reflectCoefs(i:end) = reflectCoefs(i:end) / (1 + reflectCoefs(i-1)) / (1 - reflectCoefs(i-1)) / atten;
    end
    
    gamma = calcReflect(dielEpsr(i), epsrTypes);
    magGamma = abs(gamma);
    moreDense = real(gamma) ./ abs(real(gamma));
    diffs = magGamma.*moreDense - reflectCoefs(i);
    [~,loc] = min(abs(diffs));
    
    errorRelative(i) = diffs(loc) / reflectCoefs(i);
    layers(i) = dielTypes(loc);
    dielEpsr(i+1) = epsrTypes(loc);
    
    %Following could be cleaned up
    %Depth Calculations
    [alpha(i),beta(i)] = calcGamma(dielEpsr(i)*eps0,f);
    
    dts = diff(times);
    timesInLayer = [times(1) diff(times)];
    vp = [c0, 2*pi*f./beta];
    
    depths = cumsum(vp(1:end-1).*timesInLayer) / 2;
    zs = linspace(0,d,ncells);
    ex = zeros(1,length(zs));
end

figure(1)

subplot(2,2,1)
hold on
scatter(times,rawSignal(locations),'ko')
hold off

subplot(2,2,2)
hold on
scatter(times,signal(locations),'ko')
hold off

h = 2;

figure(1)
subplot(2,2,3)
i = i+1;
simPlots(depths,layers,zs,-h,h);
hold off
