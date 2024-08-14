%% Initialization
initialization

% Path parameters
pathFolder = folder + "/RW_Paths";
pathFilenameBase = pathFolder + '/Example' + string(example) + '_' ...
                 + objective + params + "_NoDepletion";
numPaths = 1000;

% Read in files
twoDimensions = [ny, nx];
obst = readFromFile(twoDimensions, precision, filenameBase + '_Obstacle');
home = readFromFile(twoDimensions, precision, filenameBase + '_HomeBase');

downsample = 8;
nx_ds = round((nx - 1)/downsample) + 1;
ny_ds = round((ny - 1)/downsample) + 1;

% This will store frequency for a cell downsample by downsample points
cellFrequency = zeros(nx_ds - 1, ny_ds - 1);
cellFrequencyMode1 = zeros(nx_ds - 1, ny_ds - 1);
cellFrequencyMode2 = zeros(nx_ds - 1, ny_ds - 1);

h = dx;

for k = 1:numPaths
    path_number = k
    pathfile = pathFilenameBase + '_Path_' + string(k) + "_Path";
    modefile = pathFilenameBase + '_Path_' + string(k) + "_Modes";
    stepsfile = pathFilenameBase + '_Path_' + string(k) + '_Steps';

    numSteps = readFromFile(1, "int", stepsfile);
    path = readFromFile([2,numSteps], precision, pathfile);
    modes = readFromFile([numSteps, 1], 'int', modefile);
    
    for n = 1:numSteps
        pointX = path(1, n);
        pointY = path(2, n);
        mode = modes(n);
        
        % Find grid location
        floorX = floor(pointX/h);
        floorY = floor(pointY/h);

        pointI = floor(floorX/downsample) + 1;
        pointJ = floor(floorY/downsample) + 1;
        
        pointI = pointI - (pointI == nx_ds);
        pointJ = pointJ - (pointJ == ny_ds);
        
        cellFrequency(pointJ, pointI) = cellFrequency(pointJ, pointI) + 1;        

        if (mode == 1)
            cellFrequencyMode1(pointJ, pointI) = cellFrequencyMode1(pointJ, pointI) + 1;
        else
            cellFrequencyMode2(pointJ, pointI) = cellFrequencyMode2(pointJ, pointI) + 1;
        end
    end      
end

% Visualization
% Compute area normalization
numPoints = zeros(nx_ds-1, ny_ds-1);
pointsInDomain = zeros(nx_ds-1, ny_ds-1);
xx = linspace(0,1,nx);
yy = linspace(0,1,ny);
for i=1:nx
    for j=1:ny
        x = xx(i);
        y = yy(j);

        % Find grid location
        floorX = floor(x/h);
        floorY = floor(y/h);

        pointI = floor(floorX/downsample) + 1;
        pointJ = floor(floorY/downsample) + 1;
        
        pointI = pointI - (pointI == nx_ds);
        pointJ = pointJ - (pointJ == ny_ds);

        numPoints(pointJ, pointI) = numPoints(pointJ, pointI) + 1;
        if obst(j,i) == 1
            pointsInDomain(pointJ, pointI) = pointsInDomain(pointJ, pointI) + 1;
        end
    end
end

cellFrequency = cellFrequency.*(numPoints./pointsInDomain);

% Now upscale it back to the original resolution
[XX, YY] = ndgrid(xx, yy);

xx_ds = xx(downsample:downsample:end);
yy_ds = yy(downsample:downsample:end);
[XX_DS, YY_DS] = ndgrid(xx_ds, yy_ds);

% Use piecewise constant interpolation to upscale
F = griddedInterpolant(XX_DS, YY_DS, cellFrequency, 'nearest');
cellFrequencyUpscaled = F(XX, YY); 

% Remove blank space
nYShift = round(yShift*(ny-1));

% Plot results
figure();
obstacleMask = obst;
obstacleMask(obstacleMask==0) = NaN;
figureData = max(0,log(cellFrequencyUpscaled)).*obstacleMask;
pcolor(xx-h/2, yy(1:(end-nYShift))-h/2, figureData(nYShift+1:end,:)); shading flat; hold on; 
contour(xx, yy(1:(end-nYShift)), home(nYShift+1:end,:), [1 1], 'LineWidth', 2, 'LineColor', 'w');
contour(xx, yy(1:(end-nYShift)), obst(nYShift+1:end,:), [1 1], 'LineColor', 'k', 'LineWidth', 3)

% Figure Settings
ax = gca;
ax.FontSize = largelabelfontsize;
ax.FontName = figurefont;
axis equal
xlim([0,1])
ylim([0,0.7])
yticks([0 0.35 0.7])
set(gca, 'YDir', 'normal');
set(gcf,'Position', [100,100,800,600])
exportgraphics(gcf, outputFilenameBase + "_Utilization.png", 'Resolution', resolutionDPI);