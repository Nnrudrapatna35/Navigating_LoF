% Plots the optimal trajectories for various prespecified spotting times
% for the real world example. Can show with filled or not filled color
% scheme.

% Initialization
initialization;

% Path parameters
pathFolder = folder + "/RW_Paths";
pathFilenameBase = pathFolder + '/Example' + string(example) + '_' ...
                 + objective + params + "_NoDepletion";

% Additional parameters
backgroundF = false;
showFood = true;
showModeTransitions = true;

nStages = 1; %readFromFile(1, "int", filenameBase + '_NStages');

% Read files
twoDimensions = [ny, nx];
F = readFromFile(twoDimensions, precision, filenameBase + "_FoodAccumulated");
obst = readFromFile(twoDimensions, precision, filenameBase + '_Obstacle');
home = readFromFile(twoDimensions, precision, filenameBase + '_HomeBase');
mus = readFromFile(twoDimensions, precision, filenameBase + '_PredatorDensity');

% Set zoomed in area
xMin = 0;
xMax = 0.9;
yMin = 0;
yMax = 0.71;

iMin = xMin/dx + 1;
iMax = xMax/dx + 1;
jMin = yMin/dx + 1;
jMax = yMax/dx + 1;

F = F(jMin:jMax, iMin:iMax);
obst = obst(jMin:jMax, iMin:iMax);
home = home(jMin:jMax, iMin:iMax);
mus = mus(jMin:jMax, iMin:iMax);

obstacleMask = obst;
obstacleMask(obstacleMask==0) = NaN;
FMax = max(F, [], 'all');
FMin = min(F, [], 'all');
musMax = max(mus, [], 'all');
musMin = min(mus, [], 'all');

nYShift = round(yShift*(ny-1));
xx = linspace(xMin, xMax, iMax-iMin+1);
yy = linspace(yMin, yMax, jMax-jMin+1);

% Select the mode switches we would like to include
paramsList = ["", ...
              "_SpotTime_0", ...
              "_SpotTime_5", ...
              "_SpotTime_10"];
colorList = [[1 1 1];
             [1 1 1]; 
             [1 1 1]];
objectiveList = ["Linear", ...
                 "SquareRoot", ...
                 "Sigmoid"];
pathColor = colorList(objectiveList==objective, :);
lineWidth = 4.5;

% Read in trajectory files
handles = zeros(length(paramsList));
for k = 1:length(paramsList)
    % Plot the background
    figure()
    if backgroundF
        levels = linspace(FMin, FMax, 11);
        FObst = F.*obstacleMask;
        contourf(xx, yy(1:(end-nYShift)), FObst((nYShift+1):end,:), levels); hold on;
        contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
                'LineWidth', 2, 'LineColor', 'w');
    else
        levels = linspace(musMin, musMax, 9);
        musObst = mus.*obstacleMask;
        FObst = F.*obstacleMask;
    
        contourf(xx, yy(1:(end-nYShift)), musObst((nYShift+1):end,:), levels); hold on; 
        clim([musMin, musMax]);
        contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
                'LineWidth', 2, 'LineColor', 'w');
    end
    contour(xx, yy(1:(end-nYShift)), obst((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');

    % Format figure
    ax = gca;
    ax.FontSize = largelabelfontsize;
    ax.FontName = figurefont;

    
    % Read in the path
    currentParams = "_9pt_NoDepletion_NoModeSwitches" + paramsList(k);
    pathFilenameBase = pathFolder + '/Example' + string(example) ...
                     + '_' + objective + currentParams;
    
    pathfile = pathFilenameBase + "_Path";
    modefile = pathFilenameBase + "_Modes";
    energyfile = pathFilenameBase + "_Energy";
    stepsfile = pathFilenameBase + "_Steps";
    timefactorfile = pathFilenameBase + "_TimeFactor";
    
    nSteps = readFromFile(1, "int", stepsfile);
    timeFactor = readFromFile(1, "int", timefactorfile);
    maxPathLength = (timeFactor*nt - 1)*nStages + 1;
    
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
                 'Color', pathColor, 'DisplayName', objectiveName);
        else 
             plot(x(switchIndices(interval,1):switchIndices(interval,2)), ...
                 y(switchIndices(interval,1):switchIndices(interval,2))-yShift,...
                 'LineStyle', ":", 'LineWidth', lineWidth,...
                 'Color', 'k', 'HandleVisibility', 'off');
        end
    end
    
    if k~= 2
        scatter(x(1), y(1)-yShift, 125, "o", 'filled', 'MarkerFaceColor',...
                'w', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
    
    for n=2:nSteps
        if (modes(n) == 1) && (modes(n-1) == 2)
            scatter(x(n), y(n)-yShift, 150, 'o', 'filled', 'MarkerFaceColor',...
                    WongBlueGreen, 'MarkerEdgeColor', 'k', ...
                    'HandleVisibility', 'off')
        elseif (modes(n) == 2) && (modes(n-1) == 1)
                scatter(x(n), y(n)-yShift, 275, 'x', 'MarkerEdgeColor', ...
                         'k', 'LineWidth', 4.5, 'HandleVisibility', 'off')
                scatter(x(n), y(n)-yShift, 250, 'x', 'MarkerEdgeColor', ...
                        WongYellow, 'LineWidth', 4.0, 'HandleVisibility', 'off')
        end
    end
    
    if nSteps < maxPathLength
        scatter(x(nSteps), y(nSteps)-yShift, 200, 'x', 'MarkerEdgeColor',...
                WongVermillion, 'LineWidth', 4.0, 'HandleVisibility', 'off')
    end
    
    axis equal
    xlim([xMin, xMax]);
    ylim([yMin, yMax-yShift]);
    set(gcf,'Position', [100,100,1000,500])
    if backgroundF
        exportgraphics(gcf, outputFilenameBase + paramsList(k) + "_bgF_PathWindow.png");
    else
        exportgraphics(gcf, outputFilenameBase + paramsList(k) + "_PathWindow.png", 'Resolution', resolutionDPI);
    end
end

%% Plots the optimal trajectories for various prespecified spotting times
% for the real world example. Can show with filled or not filled color
% scheme.

% Initialization
initialization;
pathFolder = folder + "/RW_Paths";
pathFilenameBase = pathFolder + '/Example' + string(example) + '_' ...
                 + objective + params + "_NoDepletion";

% Additional parameters
backgroundF = false;

nStages = readFromFile(1, "int", filenameBase + '_NStages');

% Read files
twoDimensions = [ny, nx];
F = readFromFile(twoDimensions, precision, filenameBase + "_FoodAccumulated");
obst = readFromFile(twoDimensions, precision, filenameBase + '_Obstacle');
home = readFromFile(twoDimensions, precision, filenameBase + '_HomeBase');
mus = readFromFile(twoDimensions, precision, filenameBase + '_PredatorDensity');

% Set zoomed in area
xMin = 0;
xMax = 0.9;
yMin = 0;
yMax = 0.71;

iMin = xMin/dx + 1;
iMax = xMax/dx + 1;
jMin = yMin/dx + 1;
jMax = yMax/dx + 1;

F = F(jMin:jMax, iMin:iMax);
obst = obst(jMin:jMax, iMin:iMax);
home = home(jMin:jMax, iMin:iMax);
mus = mus(jMin:jMax, iMin:iMax);

obstacleMask = obst;
obstacleMask(obstacleMask==0) = NaN;
FMax = max(F, [], 'all');
FMin = min(F, [], 'all');
musMax = max(mus, [], 'all');
musMin = min(mus, [], 'all');

nYShift = round(yShift*(ny-1));
xx = linspace(xMin, xMax, iMax-iMin+1);
yy = linspace(yMin, yMax, jMax-jMin+1);

% Select the mode switches we would like to include
objectiveList = ["Linear", ...
                 "SquareRoot", ...
                 "Sigmoid"];
colorList = [[1 0 0];
             [1 1 0];
             [0 1 0]];
symbolList = ["-", ":", "-."];
lineWidth = 4.5;

% Plot the background
figure()
if backgroundF
    levels = linspace(FMin, FMax, 11);
    FObst = F.*obstacleMask;
    contourf(xx, yy(1:(end-nYShift)), FObst((nYShift+1):end,:), levels, 'HandleVisibility', 'off'); hold on;
    contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w', 'HandleVisibility', 'off');
else
    levels = linspace(musMin, musMax, 9);
    musObst = mus.*obstacleMask;
    FObst = F.*obstacleMask;

    contourf(xx, yy(1:(end-nYShift)), musObst((nYShift+1):end,:), ...
                 levels, 'HandleVisibility', 'off'); hold on; 
    clim([musMin, musMax]);
    contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w', 'HandleVisibility', 'off');
end
contour(xx, yy(1:(end-nYShift)), obst((nYShift+1):end,:), [1 1], ...
        'LineWidth', 2, 'LineColor', 'k', 'HandleVisibility', 'off');

% Read in trajectory files
handles = zeros(length(objectiveList));
for k = 1:length(objectiveList)
    
    % Read in the path
    currentParams =  "9pt_NoDepletion_NoModeSwitches";
    pathFilenameBase = pathFolder + '/Example' + string(example) ...
                     + '_' + objectiveList(k) +"_" + currentParams;

    pathfile = pathFilenameBase + "_Path";
    stepsfile = pathFilenameBase + "_Steps";

    nSteps = readFromFile(1, "int", stepsfile);

    path = readFromFile([2,nSteps], precision, pathfile);

    x = path(1,:);
    y = path(2,:);

    plot(x, y-yShift,'LineStyle', symbolList(k), 'LineWidth', lineWidth,...
         'Color', colorList(k,:), 'DisplayName', objectiveName);
    
    % Plot starting location
    scatter(x(1), y(1)-yShift, 100, "o", 'filled', 'MarkerFaceColor',...
            'w', 'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
end

% Format figure and save
ax = gca;
ax.FontSize = largelabelfontsize;
ax.FontName = figurefont;

legend(["Risk-neutral", "Risk-averse", "Sigmoid"], Location="south")

axis equal
xlim([xMin, xMax]);
ylim([yMin, yMax-yShift]);
set(gcf,'Position', [100,100,1200,600])
if backgroundF
    exportgraphics(gcf, outputFilenameBase + "_Utilities_bgF_PathWindow.png");
else
    exportgraphics(gcf, outputFilenameBase + "_Utilities_PathWindow.png", 'Resolution', resolutionDPI);
end

    
    