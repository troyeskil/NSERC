function [vertBounds, diel] = getPlotParams(d, ncells, depths, diels, zs)
h = 2;
vertBounds = [-h, h];
loc = mapToGrid(d, ncells, depths);

labels = string(zeros(1,nlayers+1));
labels(nlayers+1) = "Efield";
dielXs = zeros(1,nlayers*2);
dielYs = zeros(1,nlayers*2);

for i = 1:nlayers
    z = [zs(loc(i)),zs(loc(i))];
    dielXs(2*i-1:2*i) = z;
    dielYs(2*i-1:2*i) = vertBounds;
    labels(i) = diels(i);
end
plot(dielXs,dielYs,'--')
end