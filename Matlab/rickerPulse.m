function [signal] = rickerPulse(t, q, fc)

T = (0.8 + (1-q)/15)/fc;

vt0 = vPulse(t, T);
vt1 = vPulse(t - T/2, T);
vt2 = vPulse(t - T, T);

signal = vt0 - (2-q) * vt1 + (1-q) * vt2;
end

function [pulse] = vPulse(t1, T)

mask = (t1 > 0) .* (t1 < T);
pulse = 1/2*(1+cos(pi*(t1-T/2)/(T/2)));
pulse = pulse .* mask;
end