%% Import Files
initialization;

% Additional parameters
backgroundF = false;
showModeTransitions = true;
nStages = readFromFile(1, "int", filenameBase + '_NStages');

% Path lengths
pathFolder = folder + "/RW_Paths";
pathFilenameBase = pathFolder + '/Example' + string(example) + '_' ...
                 + objective + params + "_NoDepletion";
timeFactor = readFromFile(1, "int", pathFilenameBase + "_Path_1_TimeFactor");
maxPathLength = (timeFactor*nt-1)*nStages+1;

% Read files
twoDimensions = [ny, nx];
psiInit = readFromFile(twoDimensions, precision, initialfoodfile);
F = readFromFile(twoDimensions, precision, filenameBase + "_FoodAccumulated");
obst = readFromFile(twoDimensions, precision, filenameBase + '_Obstacle');
obstacleMask = obst;
obstacleMask(obstacleMask==0) = NaN;
home = readFromFile(twoDimensions, precision, filenameBase + '_HomeBase');
mus = readFromFile(twoDimensions, precision, filenameBase + '_PredatorDensity');

FMax = max(F, [], 'all');
FMin = min(F, [], 'all');
musMax = max(mus, [], 'all');
musMin = min(mus, [], 'all');

% Read in trajectory files
numPaths = 1000;

xx = linspace(minXY,maxXY,nx);
yy = linspace(minXY,maxXY,ny);

modes = zeros(maxPathLength, numPaths);
x = zeros(maxPathLength, numPaths);
y = zeros(maxPathLength, numPaths);
nSteps = zeros(1, numPaths);

nYShift = round(yShift*(ny-1));

figure()
if backgroundF
    levels = linspace(FMin, FMax, 11);
    FObst = F.*obstacleMask;
    contourf(xx, yy(1:(end-nYShift)), FObst((nYShift+1):end,:), levels); hold on;
    contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
else
    levels = linspace(musMin, musMax, 11);
    musObst = mus.*obstacleMask;
    contourf(xx, yy(1:(end-nYShift)), musObst((nYShift+1):end,:), levels); 
    hold on; clim([musMin, musMax]);
    contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obst((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
end


lineWidth = 1;
spottedLocations = [];
escapedLocations = [];
deathLocations = [];
        
for k = 1:numPaths
    pathfile = pathFilenameBase + '_Path_' + string(k) + "_Path";
    modefile = pathFilenameBase + '_Path_' + string(k) + "_Modes";
    energyfile = pathFilenameBase + '_Path_' + string(k) + "_Energy";
    stepsfile = pathFilenameBase + '_Path_' + string(k) + "_Steps";

    nSteps = readFromFile(1, "int", stepsfile);
    path = readFromFile([2,nSteps], precision, pathfile);
    modes = readFromFile([nSteps, 1], 'int', modefile);
    energy = readFromFile([nSteps, 1], precision, energyfile);

    x = path(1,:);
    y = path(2,:);
    mode1 = (modes==1);
    mode2 = (modes==2);
    
    % Determine when mode switches occur
    switchIndices = [1];
    modeSwitches = 0;
    for n=2:nSteps
        if (modes(n) == 1) && (modes(n-1) == 2)
            switchIndices(modeSwitches+1,2) = n;
            switchIndices(modeSwitches+2,1) = n;
            modeSwitches = modeSwitches + 1;
        elseif (modes(n) == 2) && (modes(n-1) == 1)
            switchIndices(modeSwitches+1,2) = n;
            switchIndices(modeSwitches+2,1) = n;
            modeSwitches = modeSwitches + 1;
        end
    end
    switchIndices(modeSwitches + 1,2) = nSteps;
    
    % Plot trajectory in each mode
    for interval = 1:modeSwitches+1
        if mod(interval,2) == 1
            plot(x(switchIndices(interval,1):switchIndices(interval,2)), ...
                 y(switchIndices(interval,1):switchIndices(interval,2))-yShift,...
                 'LineStyle', "-", 'LineWidth', lineWidth,...
                 'Color', 'w', 'HandleVisibility', 'off');
            
            if (switchIndices(interval,2) ~= nSteps)
                spottedLocations = [spottedLocations; x(switchIndices(interval,2)), y(switchIndices(interval,2))-yShift];
            end
            
        else 
             plot(x(switchIndices(interval,1):switchIndices(interval,2)), ...
                 y(switchIndices(interval,1):switchIndices(interval,2))-yShift,...
                 'LineStyle', "-.", 'LineWidth', lineWidth,...
                 'Color', 'k', 'HandleVisibility', 'off');

             if (switchIndices(interval,2) ~= nSteps)
                escapedLocations = [escapedLocations; x(switchIndices(interval,2)), y(switchIndices(interval,2))-yShift];
            end
        end
    end

    % Plot mode switches
    scatter(x(1), y(1)-yShift, 75, '.', 'filled', 'MarkerFaceColor',...
            'k', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off')

    if nSteps < maxPathLength
        deathLocations = [deathLocations; x(end), y(end)-yShift];
    end
end

if ~isempty(spottedLocations)
    scatter(spottedLocations(:,1), spottedLocations(:,2), 100, 'x', ...
            'MarkerEdgeColor', WongYellow, 'LineWidth', 3.0, ...
            'HandleVisibility', 'off')
    
    scatter(escapedLocations(:,1), escapedLocations(:,2), 'o', 'filled', ...
            'MarkerFaceColor', WongBlueGreen, 'MarkerEdgeColor', WongBlueGreen, ...
            'HandleVisibility', 'off')

    scatter(deathLocations(:,1), deathLocations(:,2), 100, 'x', ...
            'MarkerEdgeColor',WongVermillion, 'LineWidth', 3.0, ...
            'HandleVisibility', 'off')
end

ax = gca;
ax.FontSize = largelabelfontsize;
ax.FontName = figurefont;

axis equal
xlim([0 1]);
ylim([0, 0.7]);
yticks([0 0.35 0.7]);
set(gcf,'Position', [100,100,800,600])

if backgroundF
    exportgraphics(gcf, outputFilenameBase + "_bgF_PathFrequency.png");
else
    exportgraphics(gcf, outputFilenameBase + "_PathFrequency.png", "Resolution", resolutionDPI);
end

    