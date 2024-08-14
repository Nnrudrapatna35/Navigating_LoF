%% Plot one path at a time
initialization;

% Additional Parameters
numPaths = 6;
pathLabel = "_realization_";
mus = 3;
muk = 3;
arrowPosition = 40;
arrowSize = 4.5;

% Figure settings
lineWidth = 3;
alphaValue = 0.9;
mode2Colors = zeros(6,3);

timeFactor = readFromFile(1, "int", filenameBase + "_Stage0_TimeFactor");
nStages = readFromFile(1, "int", filenameBase + '_NStages');
maxPathLength = (timeFactor*(nt-1)-1)*nStages+1; 
stageEndIndices = ([1:nStages]*(maxPathLength - 1)/(nStages))+1;

% Read in background files
twoDimensions = [ny, nx];
F = readFromFile(twoDimensions, precision, filenameBase + "_FoodAccumulated");
obst = readFromFile(twoDimensions, precision, filenameBase + '_Obstacle');
home = readFromFile(twoDimensions, precision, filenameBase + '_HomeBase');

obstacleMask = obst;
obstacleMask(obstacleMask==0) = NaN;
FMax = max(F, [], 'all');
FMin = min(F, [], 'all');

xx = linspace(minXY,maxXY,nx);
yy = linspace(minXY,maxXY,ny);

% Plot background data
levels = linspace(FMin, FMax, 7);
FObst = F.*obstacleMask;

for k = 1:numPaths
    spottedLocations = [];
    escapedLocations = [];

    % Plot background
    figure()
    contourf(xx, yy(1:(end-nYShift)), FObst((nYShift+1):end,:), levels); hold on;
    contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2.75, 'LineColor', 'k');
    contour(xx, yy(1:(end-nYShift)), home((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2.5, 'LineColor', 'w');

    % Read in path
    currentParams = "_9pt_N_5_mus_" + string(mus) + "_muk_" + string(muk) ...
                 + pathLabel + string(k-1);
    pathFilenameBase = folder + '/Example' + string(example)...
                        + '_' + objective + currentParams;

    pathfile = pathFilenameBase + "_OptimalTrajectory";
    modefile = pathFilenameBase + "_ModeList";
    stepsfile = pathFilenameBase + "_TotalSteps";

    nSteps = readFromFile(1, "int", stepsfile);
    path = readFromFile([2,nSteps], precision, pathfile);
    modes = readFromFile([nSteps, 1], 'int', modefile);
    nSurvivedStages = readFromFile(1, "int", pathFilenameBase + "_NSurvivedStages");
    
    % Identify indexes of mode switches
    x = path(1,:);
    y = path(2,:);
    mode1 = (modes==1);
    mode2 = (modes==2);

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
    
    % Plot each mode interval
    for interval = 1:modeSwitches+1
        if mod(interval,2) == 1
            mode1Shift = 0.000*(rand(1,1) - 0.5);
            plot(x(switchIndices(interval,1):switchIndices(interval,2)) + mode1Shift,...
                 y(switchIndices(interval,1):switchIndices(interval,2))-yShift + mode1Shift,...
                 'LineStyle', '-', 'LineWidth', (1.2 - 0.05*interval)*lineWidth, ...
                 'Color', [1 1 1 alphaValue-0.1]);
            if (switchIndices(interval,2)-switchIndices(interval,1) > timeFactor*arrowPosition+arrowSize*timeFactor)

                arrowStart = path(:, round(switchIndices(interval,1)+timeFactor*arrowPosition));
                arrowEnd = path(:, round(switchIndices(interval,1)+timeFactor*arrowPosition) + arrowSize*timeFactor);
                axis manual;
                norm(arrowStart-arrowEnd)
                if norm(arrowStart - arrowEnd) > 0.005
                    arrow(arrowStart, arrowEnd, 'Color', 'c','Length', 25);
                end
            end
        else 
            mode2Shift = 0.000*(rand(1,1) - 0.5);
            plot(x(switchIndices(interval,1):switchIndices(interval,2)) + mode2Shift, ...
                 y(switchIndices(interval,1):switchIndices(interval,2))-yShift + mode2Shift,...
                 'LineStyle', ':', 'LineWidth', (1.5-0.05*interval)*lineWidth,...
                 'Color', [mode2Colors(interval/2,:) alphaValue]);
            if (switchIndices(interval,2)-switchIndices(interval,1) > timeFactor*arrowPosition+arrowSize*timeFactor)

                arrowStart = path(:, round(switchIndices(interval,1)+timeFactor*arrowPosition));
                arrowEnd = path(:, round(switchIndices(interval,1)+timeFactor*arrowPosition) + timeFactor);
                axis manual;
                norm(arrowStart - arrowEnd)
                if norm(arrowStart - arrowEnd) > 0.0025
                    arrow(arrowStart, arrowEnd, 'Color', 'c','Length', 25);
                end
            end
        end
    end

    % Extra arrows to help visualize behavior
    if objective == "SquareRoot" || (objective == "Sigmoid" && k ~= 3) || objective == "Linear"
        if nSurvivedStages >= 4
            stage4Start = stageEndIndices(3);
            arrowStart = path(:, round(stage4Start+timeFactor*arrowPosition));
            arrowEnd = path(:, round(stage4Start+timeFactor*arrowPosition) + timeFactor);
            axis manual;
            norm(arrowStart - arrowEnd)
            if norm(arrowStart - arrowEnd) > 0.0025
                arrow(arrowStart, arrowEnd, 'Color', 'c','Length', 25);
            end
        end
    end
    if (objective == "Sigmoid" && k >= 3) || (objective == "Linear" && k == 3)
        if nSurvivedStages >= 4
            stage5Start = stageEndIndices(4);
            arrowStart = path(:, round(stage5Start+timeFactor*arrowPosition));
            arrowEnd = path(:, round(stage5Start+timeFactor*arrowPosition) + timeFactor);
            axis manual;
            norm(arrowStart - arrowEnd)
            if norm(arrowStart - arrowEnd) > 0.0025
                arrow(arrowStart, arrowEnd, 'Color', 'c','Length', 25);
            end
        end
    end
    if objective == "Sigmoid" && k == 6
        stage4End = stageEndIndices(4);
        arrowStart = path(:, round(stage4End-timeFactor*4*arrowPosition) - timeFactor);
        arrowEnd = path(:, round(stage4End-timeFactor*4*arrowPosition));
        axis manual;
        norm(arrowStart - arrowEnd)
        if norm(arrowStart - arrowEnd) > 0.0025
            arrow(arrowStart, arrowEnd, 'Color', 'c','Length', 25);
        end

        stage5End = stageEndIndices(5);
        arrowStart = path(:, round(stage5End-timeFactor*2*arrowPosition) - timeFactor);
        arrowEnd = path(:, round(stage5End-timeFactor*2*arrowPosition));
        axis manual;
        norm(arrowStart - arrowEnd)
        if norm(arrowStart - arrowEnd) > 0.0025
            arrow(arrowStart, arrowEnd, 'Color', 'c','Length', 25);
        end
    end


    % Plot death locations
    if nSteps < maxPathLength
        scatter(x(nSteps), y(nSteps)-yShift, 100, 'x', 'MarkerEdgeColor',...
                WongVermillion, 'LineWidth', 4.0, 'HandleVisibility', 'off')
    end

    grey = [0.5 0.5 0.5];
    scatter(x(1), y(1), 'o', 'filled', 'LineWidth', 1.5, ...
            'MarkerEdgeColor', grey, 'MarkerFaceColor', grey)
    currentStageEndIndices = stageEndIndices(1:nSurvivedStages);
    stageEndColor = 0.8*WongRedPurple+0.2*[1 0 1];
    scatter(x(currentStageEndIndices), y(currentStageEndIndices), 'o', 'filled', 'LineWidth', 1, ...
            'MarkerEdgeColor', stageEndColor, 'MarkerFaceColor', stageEndColor);
    
    ax = gca;
    ax.FontSize = labelfontsize;
    ax.FontName = figurefont;
    
    axis equal
    xlim([minXY, maxXY]);
    ylim([minXY, maxXY-yShift]);
    yticks([0 0.2 0.4 0.6])
    set(gcf,'Position', [100,100,400,300])
    outputFilename = "../plots/" + outputFolder + '/Example' + string(example)...
                     + "_N_5_mus_" + string(mus) + "_muk_" + string(muk) ...
                     + "_" + objective + pathLabel + string(k-1);
    
    exportgraphics(gcf, outputFilename+ "_PathImage.png", 'Resolution', resolutionDPI);
end

%% Plot different visualization
initialization;

% Additional Parameters
numPaths = 6;
pathLabel = "_realization_";
mus = 3;
muk = 3;
arrowPosition = 45;
arrowSize = 5;

% Figure settings
lineWidth = 3;
alphaValue = 0.9;
mode2Colors = zeros(6,3);

timeFactor = readFromFile(1, "int", filenameBase + "_Stage0_TimeFactor");
nStages = readFromFile(1, "int", filenameBase + "_NStages")
maxPathLength = (timeFactor*(nt-1)-1)*nStages+1; 
stageEndIndices = ([1:nStages]*(maxPathLength - 1)/(nStages))+1;

% Read in background files
twoDimensions = [ny, nx];
psiInit = readFromFile(twoDimensions, precision, initialfoodfile);
F = readFromFile(twoDimensions, precision, filenameBase + "_FoodAccumulated");
obst = readFromFile(twoDimensions, precision, filenameBase + '_Obstacle');
home = readFromFile(twoDimensions, precision, filenameBase + '_HomeBase');

obstacleMask = obst;
obstacleMask(obstacleMask==0) = NaN;
FMax = max(F, [], 'all');
FMin = min(F, [], 'all');

xx = linspace(minXY,maxXY,nx);
yy = linspace(minXY,maxXY,ny);

% Plot background data
levels = linspace(FMin, FMax, 7);
FObst = F.*obstacleMask;
grey = [0.5 0.5 0.5];


for k = 1:numPaths
    figure()
    hold on;
    
    distanceHomeCorner = sqrt(2)*0.1/2;
    distanceFoodCenter = sqrt((0.4-0.1)^2 + (0.3-0.1)^2);

    % Read in path
    currentParams = "_9pt_N_5_mus_" + string(mus) + "_muk_" + string(muk) ...
                 + pathLabel + string(k-1);
    pathFilenameBase = folder + '/Example' + string(example)...
                        + '_' + objective + currentParams;

    pathfile = pathFilenameBase + "_OptimalTrajectory";
    modefile = pathFilenameBase + "_ModeList";
    stepsfile = pathFilenameBase + "_TotalSteps";

    nSteps = readFromFile(1, "int", stepsfile);
    path = readFromFile([2,nSteps], precision, pathfile);
    modes = readFromFile([nSteps, 1], 'int', modefile);
    nSurvivedStages = readFromFile(1, "int", pathFilenameBase + "_NSurvivedStages");

    % Plot home and food distances
    t = linspace(0, nSteps/maxPathLength, nSteps);
    plot(t, distanceHomeCorner*ones(size(t)), 'LineWidth', 0.5*lineWidth, ...
             'Color', grey)
    plot(t, distanceFoodCenter*ones(size(t)), 'LineWidth', 0.5*lineWidth, ...
             'Color', WongBlueGreen)
    
    % Identify indexes of mode switches
    x = path(1,:);
    y = path(2,:);
    mode1 = (modes==1);
    mode2 = (modes==2);

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

    signedDistance = sqrt((0.1 - path(1,:)).^2 + (0.3 - path(2,:)).^2);
    pointLabels = sign(path(2,:) - 0.3).*modes';
    changePoints = [0 find(diff(pointLabels)) nSteps-1] + 1;

    for index=2:length(changePoints)
        sectionStart = changePoints(index-1);
        sectionEnd = changePoints(index);
        lastMode=pointLabels(sectionStart);
        if lastMode == 1
            pathColor = (WongVermillion + [0.8 0.1 0.3])/2;
            pathStyle = '-.';
        elseif lastMode == -1
            pathColor = WongOrange;
            pathStyle = '-';
        else
            pathColor = 'k';
            pathStyle = ':';
        end

        plot(t(sectionStart:sectionEnd), signedDistance(sectionStart:sectionEnd),...
            'LineStyle', pathStyle, 'LineWidth', 0.9*lineWidth, ...
             'Color', pathColor);
    end

    % Plot death locations
    if nSteps < maxPathLength
        scatter(t(nSteps), signedDistance(nSteps)-yShift, 100, 'x', 'MarkerEdgeColor',...
                WongVermillion, 'LineWidth', 4.0, 'HandleVisibility', 'off')
    end

    scatter(t(1), signedDistance(1), 'o', 'filled', 'LineWidth', 1.5, ...
            'MarkerEdgeColor', grey, 'MarkerFaceColor', grey)
    currentStageEndIndices = stageEndIndices(1:nSurvivedStages);
    stageEndColor = 0.8*WongRedPurple+0.2*[1 0 1];
    scatter(t(currentStageEndIndices), signedDistance(currentStageEndIndices), 'o', 'filled', 'LineWidth', 1, ...
            'MarkerEdgeColor', stageEndColor, 'MarkerFaceColor', stageEndColor)
    
    ax = gca;
    ax.FontSize = labelfontsize;
    ax.FontName = figurefont;
    
    axis equal
    xlim([0, 1]);
    ylim([0, 0.4]);
    set(gcf,'Position', [100,100,750,300])
    outputFilename = "../plots/" + outputFolder + '/Example' + string(example)...
                     + "_N_5_mus_" + string(mus) + "_muk_" + string(muk) ...
                     + "_" + objective + pathLabel + string(k-1);
    
    exportgraphics(gcf, outputFilename+ "_PathOverTime.png", 'Resolution', resolutionDPI);
end