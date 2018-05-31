function [gamma] = calcReflect(epsr1, epsr2)
gamma = (sqrt(epsr1) - sqrt(epsr2)) ./ (sqrt(epsr1) + sqrt(epsr2));