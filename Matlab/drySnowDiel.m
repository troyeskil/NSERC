function [eps] = drySnowDiel(snowDensity, Tk, f)

epsIce = pureIceDiel(Tk, f);
vi = snowDensity / 0.9167;
eps = 1 + 3*vi*(epsIce - 1) / ((2+epsIce) - vi*(epsIce-1));
end