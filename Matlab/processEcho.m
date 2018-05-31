%function [diels, depths] = processEcho(t, rawSignal, f, sentPulse)

%close all
clear all
clc

i = 1;
i = 1;
powerThreshold = 0.005;
f = 1e9;

%imports
rawSignal = csvread('echoSignal.csv');
t = csvread('timeStamps.csv');
sentPulse = csvread('sentPulse.csv');

figure(1)
subplot(4,4,i)
i = i+1;
plot(t,rawSignal);

%Initial Calculations
pulseT = 1/f;
nsteps = length(t);


%check cross correlation
r = xcorr(rawSignal,sentPulse);
%r = r(1:nsteps);

% figure(1)
% subplot(4,4,i)
% i = i+1;
% plot(r);
% title('Xcorr')
% 
% deriv = diff(sentPulse);
% 
% sentAvg = nPointAvg(sentPulse, 7);
% derivAvg = nPointAvg(deriv,7);
% 
% subplot(4,4,i)
% i = i+1;
% plot(sentPulse);
% title('Sent')
% 
% subplot(4,4,i)
% i = i+1;
% plot(deriv);
% title('Sent deriv')
% 
% subplot(4,4,i)
% i = i+1;
% plot(sentAvg);
% title('Sent Avg')
% 
% subplot(4,4,i)
% i = i+1;
% plot(derivAvg);
% title('Sent deriv Avg')
% 
 corrLength = length(sentPulse) + length(rawSignal) - 1;
% PulseF = fft(sentPulse,corrLength);
% RawF = fft(rawSignal,corrLength);
% 
% subplot(4,4,i)
% i = i+1;
% plot(fftshift(abs(PulseF)));
% title('PulseF')
% 
% subplot(4,4,i)
% i = i+1;
% plot(fftshift(abs(RawF)));
% title('RawF')
% 
% Hf = conj(PulseF).*RawF;
% 
% r = fftshift(ifft(Hf));
% 
% yt = r - xcorr(sentPulse, sentPulse);
% Yf = fft(yt);
% 
% Xf = Yf ./ conj(PulseF);
% xt = ifft(Xf);
% 
% Hprime = xcorr(xt,sentPulse);
% 
% Mf = fft(r);
% Nf = Mf ./ conj(PulseF);
% nt = ifft(Nf);
% 
% subplot(4,4,i)
% i = i+1;
% plot(yt);
% title('yt')
% 
% subplot(4,4,i)
% i = i+1;
% plot(fftshift(abs(Yf)));
% title('Yf')
% 
% subplot(4,4,i)
% i = i+1;
% plot(fftshift(abs(Xf)));
% title('Xf')
% 
% subplot(4,4,i)
% i = i+1;
% plot(t,xt(222:end));
% title('xt')
% 
% subplot(4,4,i)
% i = i+1;
% plot(fftshift(abs(Mf)));
% title('Mf')
% 
% subplot(4,4,i)
% i = i+1;
% plot(fftshift(abs(Nf)));
% title('Nf')
% 
% subplot(4,4,i)
% i = i+1;
% plot(t,nt(222:end));
% title('nt')
% 
% subplot(4,4,i)
% i = i+1;
% plot(Hprime(443:end));
% title("H't")
% 
% subplot(4,4,i)
% i = i+1;
% plot(inverseXcorr(sentPulse,r));
% title("a")

%Processing of Signal
locPulse = findExtrema(1,nsteps,sentPulse);
shiftPulse = locPulse - 1;

signal = rawSignal - sentPulse;
r = xcorr(signal, sentPulse);
keepLoop = 1;
echoPower = sum(signal.^2);
thresh = echoPower / 100;
j = 1;
k = 2;
ratio(1) = 1;
timeDelay(1) = 0;
aboveThresh = zeros(1,corrLength);
%signal(nsteps) = 1;


%using autocorrelation for echo cancellation
% f5 = 1;
% minPk = 0.05;
% 
% [Rmm, lags] = xcorr(rawSignal, 'unbiased');
% 
% Rmm2 = Rmm(lags>0);
% 
% [pks,d1] = findpeaks(Rmm2,lags(lags>0),'MinPeakHeight', 0.05);
% 
% figure(5)
% subplot(4,4,f5)
% f5 = f5+1;
% plot(Rmm)
% title("Rmm")
% 
% subplot(4,4,f5)
% f5 = f5+1;
% plot(lags)
% title("lags")
% 
% subplot(4,4,f5)
% f5 = f5+1;
% plot(Rmm2)
% title("Rmm2")
% 
% subplot(4,4,f5)
% f5 = f5+1;
% plot(lags>0)
% title("lags")
minPk = 0.5;

while keepLoop
    oldSignal = signal;
    signalEnergy = signal.^2;
    aboveThresh = abs(signal) > thresh;
    loc(1) = find(aboveThresh,1);
    loc(2) = loc(1) - 1 + find(~aboveThresh(loc(1):end),1);
    
    if (k > 4 || loc(1) == nsteps || loc(1) > 0.9*nsteps)
        keepLoop = 0;
    else
        echoLoc = findExtrema(loc(1),loc(2),signal);
        Emag = signal(echoLoc);
        ratio(k) = Emag / sentPulse(locPulse);
        timeDelay(k) = t(echoLoc - shiftPulse);

        %%%%%
        %Time domain echo Removal
        %%%%%
        echo = zeros(1,nsteps);
        shift = echoLoc - locPulse;
        lengthLeft = nsteps - (echoLoc - shiftPulse);
        echo(echoLoc - shiftPulse: end) = sentPulse(1:lengthLeft + 1) * ratio(k);
        signal = signal - echo;
        
        %%%%%
        %Echo removal using cross correlation
        %%%%%
        [Rmm, lags] = xcorr(oldSignal, sentPulse);
        lengthcorr = length(Rmm);
        Rmm2 = Rmm(lags>0);
        posLags = lags(lags>0);
        [absPks,locs] = findpeaks(abs(Rmm2),'MinPeakHeight', minPk);
        pks = Rmm2(locs);
        
        [rPulse, ~] = xcorr(sentPulse, sentPulse);
        rPulse2 = rPulse(lags>0);
        [pulseAbsPks, pulseLocs] = findpeaks(abs(rPulse2), 'MinPeakHeight', minPk);
        pulsePks = rPulse2(pulseLocs);
        
        mid = find(lags>0,1);
        scale = pks(1) / pulsePks(end);
        shift = locs(1) - (mid - pulseLocs(end));
        shift = -shift;
        
        rEcho = zeros(1,mid-2);
        rEcho = rPulse(shift:shift + mid - 3) * scale;
        
        newRmm2 = Rmm2 - rEcho;
        newRmm(1:mid-1) = 0;
        newRmm(mid:lengthcorr) = newRmm2;
        
        newSignal = inverseXcorr(sentPulse,newRmm);
        old = inverseXcorr(sentPulse,Rmm);
        

        keepLoop = 0;
        k = k+1;
        j = 1;
        
        %Plot
        figure(k-1)

        subplot(4,4,j)
        j = j+1;
        plot(oldSignal)
        title("Signal")
        
        subplot(4,4,j)
        j = j+1;
        plot(signal)
        title("Signal w/ echo Removed")
        
        subplot(4,4,j)
        j = j+1;
        plot(oldSignal.^2)
        title("Signal Power")

        subplot(4,4,j)
        j = j+1;
        plot(echo)
        title("Echo")
        
        subplot(4,4,j)
        j = j+1;
        plot(Rmm)
        title("Rmm")
        
        subplot(4,4,j)
        j = j+1;
        plot(Rmm2)
        title("Rmm2")
        hold on
        scatter(locs,pks,'ko')
        hold off
        
        subplot(4,4,j)
        j = j+1;
        plot(rPulse)
        title("rPulse")
        
        subplot(4,4,j)
        j = j+1;
        plot(rPulse2)
        title("rPulse2")
        hold on
        scatter(pulseLocs,pulsePks,'ko')
        hold off
        
        subplot(4,4,j)
        j = j+1;
        plot(rEcho)
        title("rEcho")
        
        subplot(4,4,j)
        j = j+1;
        plot(newRmm2)
        title("newRmm2")
        
        subplot(4,4,j)
        j = j+1;
        plot(newRmm)
        title("newRmm")
        
        subplot(4,4,j)
        j = j+1;
        plot(newSignal)
        title("newSignal")
        
        subplot(4,4,j)
        j = j+1;
        plot(old)
        title("Rmm inverse xcorr")
        
        
        subplot(4,4,j)
        j = j+1;
        plot(aboveThresh)
        title("Above Threshold")
        
        thresholdLine = ones(1,corrLength)*thresh;
        subplot(4,4,j)
        j = j+1;
        plot(oldSignal)
        hold on
        plot(thresholdLine)
        plot(-thresholdLine)
        hold off
        
    end
end

sentPower = sum(sentPulse.^2);
nrcs = echoPower / sentPower;

figure(1)
thresholdLine = ones(1,corrLength)*thresh;
subplot(4,4,i)
i = i+1;
plot(r)
hold on
plot(thresholdLine)
hold on
plot(-thresholdLine)
hold off

for k = 2:length(ratio)
    
end

