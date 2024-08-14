% Read in data as a matrix and interpolate/downsample to get a data matrix
% of desired dimensions

data = "MuS";
example = "7";
filename = data + ".csv";

grid = readmatrix(filename);
oldSize = size(grid);
nxOld = oldSize(1);
nyOld = oldSize(2);
% Test region size 1
% xx = linspace(0, 1.4, nxOld);
% yy = linspace(0, 1.4, nyOld);
% Test region size 2
physMax = 2;
xx = linspace(0, physMax, nxOld);
yy = linspace(0, physMax, nyOld);

[XX, YY] = ndgrid(xx, yy);

F = griddedInterpolant(XX, YY, grid);

nxNew = 201;
nyNew = 201;

% Test region 1
% xxNew = linspace(0.4, 1.4, nxNew);
% yyNew = linspace(0, 0.95, nyNew);
% Whole region
xxNew = linspace(0, physMax, nxNew);
yyNew = linspace(0, physMax, nyNew);

[XXNew, YYNew] = ndgrid(xxNew, yyNew);
newData = F(XXNew, YYNew);

if strcmp(data, "Obstacle") == 1
    newData = floor(newData); % Round conservatively to determine what is obstacle
elseif strcmp(data, "SleepingsSites") == 1
    newData = floor(newData); % Round conservatively to determine what is obstacle
end

contourf(grid)
figure
contourf(newData)


mkdir("Example" + example)
outputFilename = "Example" + example + "/" + data + "_" + string(nxNew) + ".csv";

writematrix(newData, outputFilename)