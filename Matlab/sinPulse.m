function [signal] = sinPulse(t, f)
    period = 1/f;
    maskSinPulse = t < period;
    signal = sin(2*pi*f*t) .* maskSinPulse;
end