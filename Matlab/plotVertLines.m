function [] = plotVertLines(xs,locs,hiBound,loBound)
n = length(xs);
bound = [loBound,hiBound];
for loc = locs
    i = mapToGrid(xs(end),n,loc);
    x = [xs(i) xs(i)];
    plot(
    