close all
clear all
clc

f = 1e9;
saveRCS = 1;
saveSignal = 1;
viewMovie = 0;
saveMovie = 1;
speedFactor = 1;
fps = 12;

for scenario = 1:5
    FDTDv5(scenario,f,saveRCS,saveSignal,viewMovie,saveMovie,speedFactor,fps)
    close all
end