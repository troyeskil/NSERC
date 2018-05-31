function [i] = mapToGrid(gridLength, gridPoints, loc)

pointsPerMeter = gridPoints / gridLength;
i = ceil(pointsPerMeter * loc) + 1;
end