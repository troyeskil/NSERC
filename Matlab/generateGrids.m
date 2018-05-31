function [epsRe, condEff, epsImag, diels] = generateGrids(d, ncells, depths, diels, epsrSpecified, Tk, f, snowDensity)
eps0 = 8.854e-12;
epsr = ones(1,ncells);
S = 0.3254;
for i = 1:length(depths)
    layer = diels(i);
    if (layer == "snow")
        epsLayer = drySnowDiel(snowDensity, Tk(i), f);
    elseif (layer == "ice")
        epsLayer = pureIceDiel(Tk(i), f);
    elseif (layer == "air")
        epsLayer = 1;
    elseif (layer == "sea water")
        T = Tk(i);
        if Tk(i) <= 273.15
            T = 273.15;
        end
        epsLayer = seaWaterDiel(f,T,S);
    elseif (layer == "water")
        T = Tk(i);
        if Tk(i) <= 273.15
            T = 273.15;
        end        
        epsLayer = waterDiel(f,T);
    elseif (layer == "brine")
        epsLayer = brineDiel(f,Tk(i));
    else
        epsLayer = epsrSpecified(i);
        diels(i) = "Dielectric with epsr = " + num2str(epsLayer);
    end
    startCell = mapToGrid(d,ncells, depths(i));
    epsr(startCell:ncells) = epsLayer;
end
epsImag = imag(epsr) * eps0;
epsRe = real(epsr) * eps0;
condEff = 2 * pi * f * imag(epsr) * eps0;
end