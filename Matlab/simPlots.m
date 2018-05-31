function [p] = simPlots(depths, diels, zs, loBound, hiBound)

vertBounds = [loBound, hiBound];
d = zs(end);
ncells = length(zs);
nlayers = length(depths);

ex = zeros(1,ncells);
loc = mapToGrid(d, ncells, depths);

labels = string(zeros(1,nlayers+1));
labels(1) = "Efield";

p = plot(zs, ex);
ylim(vertBounds);
hold on

for i = 1:nlayers
    z = [zs(loc(i)),zs(loc(i))];
    plot(z, vertBounds,'--')
    hold on
    labels(i+1) = diels(i);
end

legend(labels)
hold on
end