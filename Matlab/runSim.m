close all
clear all
clc

f = 1e9;
saveRCS = 0;
saveSignal = 1;
viewMovie = 0;
saveMovie = 0;
speedFactor = 10;
fps = 24;
scenario = 6;
attenuateLength = 1;
ratio = 1.6/5;

recieved = FDTDv5(scenario,f,saveRCS,saveSignal,viewMovie,saveMovie,speedFactor,fps);

if attenuateLength
    rawSignal = csvread('echoSignal.csv');
    t = csvread('timeStamps.csv');
    sentPulse = csvread('sentPulse.csv');
    
    length = ceil(length(rawSignal) * ratio);
    echo = rawSignal(1:length);
    t2 = t(1:length);
    Esource = sentPulse(1:length);
    
    csvwrite('echoSignal.csv',echo);
    csvwrite('timeStamps.csv',t2);
    csvwrite('sentPulse.csv', Esource);
end
