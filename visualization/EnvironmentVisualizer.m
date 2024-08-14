% Generate contour plots of environment characteristics: metabolic cost,
% speed, mode switching rates, and food density/availability

initialization;

% Global figure settings
filled_color = true;

% Set the number of level sets based on the example
if example == "RealWorld"
    nlevels = 11;
    labelfontsize = largelabelfontsize;
else
    nlevels = 7;
end

% Read in background information
xx = linspace(minXY, maxXY, nx);
yy = linspace(minXY, maxXY, ny);

obstacle = readFromFile([ny nx], precision, filenameBase + '_Obstacle');
obstacleMask = obstacle;
obstacleMask(obstacleMask==0) = NaN;
homeBase = readFromFile([ny nx], precision, filenameBase + '_HomeBase');

%% Plot metabolic cost
datafilename = filenameBase + '_Stage' + string(stageNum) + '_Mode' ...
               + string(modeNum) + 'Cost';
K = readFromFile([ny nx], precision, datafilename);
Kmin = min(K.*obstacle, [], 'all');
Kmax = max(K.*obstacle, [], 'all');

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), K((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
else
    contourf(xx, yy, K.*obstacle, nlevels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2, 'LineColor', 'w');
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;
cbr = colorbar;
if Kmin == Kmax
    cbr.TicksMode = 'Manual';
    cbr.Ticks = Kmin;
    tix = cbr.Ticks;
    cbr.TickLabels = compose('%9.2f',tix);   
end

axis equal
plotfilename = strcat(outputFilenameBase + '_Mode' + string(modeNum) ...
               + '_Cost.png');
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end
exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);

%% Plot Speed
datafilename = filenameBase + '_Stage' + string(stageNum) + '_Mode' ...
             + string(modeNum) + 'Speed';
f = readFromFile([ny nx], precision, datafilename);

fmin = min(f.*obstacle, [], 'all');
fmax = max(f.*obstacle, [], 'all');
levels = linspace(fmin, fmax, nlevels);

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), f((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
else
    contourf(xx, yy, f.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2, 'LineColor', 'w');
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;
cbr = colorbar;
if fmin == fmax
    cbr.TicksMode = 'Manual';
    cbr.Ticks = fmin;
    tix = cbr.Ticks; 
    cbr.TickLabels = compose('%9.2f',tix); 
end
 

plotfilename = strcat(outputFilenameBase + '_Mode' + string(modeNum) ...
             + '_Speed.png');
axis equal
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end
exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);

%% Plot Spotting Rate
datafilename = filenameBase + '_PredatorDensity';
muS = readFromFile([ny nx], precision, datafilename);
musmin = min(muS.*obstacle, [], 'all');
musmax = max(muS.*obstacle, [], 'all');

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), muS((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
elseif example == "MultiStage"
    muS = muS.*obstacle;
    musmin = min(muS(homeBase == 1), [], 'all');
    levels = [0 linspace(musmin, musmax, 5)];
    contourf(xx, yy, muS.*obstacle, levels); hold on; 
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.75, 'LineColor', 'k');
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.5, 'LineColor', 'w');
    yticks([0 0.2 0.4 0.6])
    clim([musmin/2 musmax])
else
    if filled_color
        contourf(xx, yy, muS.*obstacle, nlevels); hold on;
        contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.2, 'LineColor', 'w');
    else
        contour(xx, yy, muS.*obstacle, nlevels, 'LineWidth', 2); hold on;
        colormap("autumn")
        contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.2, 'LineColor', WongBlue);
        clim([musmin musmax])
    end
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;
cbr = colorbar;
if musmin == musmax
    cbr.TicksMode = 'Manual';
    cbr.Ticks = musmin;
end
cbr.Ticks
cbr.TickLabels = compose('%.2f',cbr.Ticks); 

plotfilename = strcat(outputFilenameBase + '_SpottingRate.png');
axis equal
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end
cbr.Ticks
cbr.TickLabels = compose('%.2f',cbr.Ticks); 

exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);

%% Plot kill rate
datafilename = filenameBase + '_KillRate';
muK = readFromFile([ny nx], precision, datafilename);
mukmin = min(muK.*obstacle, [], 'all');
mukmax = max(muK.*obstacle, [], 'all');

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), muK((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
elseif example == "MultiStage"
    muK = muK.*obstacle;
    mukmin = min(muK(homeBase == 1), [], 'all');
    levels = [0 linspace(mukmin, mukmax, 5)];
    contourf(xx, yy, muK.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.75, 'LineColor', 'k');
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.5, 'LineColor', 'w');
    yticks([0 0.2 0.4 0.6])
    clim([mukmin/2 mukmax])
else
    contourf(xx, yy, muK.*obstacle, nlevels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2, 'LineColor', 'w');
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;

plotfilename = strcat(outputFilenameBase + '_KillRate.png');
axis equal
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end

cbr = colorbar;
if mukmin == mukmax
    cbr.TicksMode = 'Manual';
    cbr.Ticks = mukmin;
end
cbr.TickLabels = compose('%.1f',cbr.Ticks);
exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);

%% Plot give up rate
datafilename = filenameBase + '_GiveUpRate';
muG = readFromFile([ny nx], precision, datafilename);

mugmin = min(muG.*obstacle, [], 'all');
mugmax = max(muG.*obstacle.*(homeBase), [], 'all');

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), muG((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
elseif example == "MultiStage"
    levels = [mugmin-0.5 linspace(mugmin, mugmax, 7) mugmax+0.5];
    contourf(xx, yy, muG.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.75, 'LineColor', 'k');
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.5, 'LineColor', 'w');
    clim([mugmin-0.1 mugmax+0.1])
    yticks([0 0.2 0.4 0.6])
else
    levels = linspace(mugmin, mugmax, nlevels);
    contourf(xx, yy, muG.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2, 'LineColor', 'w');
    clim([mugmin mugmax])
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;
cbr = colorbar;
% clim([mugmin mugmax])
if mugmin == mugmax
    cbr.TicksMode = 'Manual';
    cbr.Ticks = mugmin;
end

plotfilename = strcat(outputFilenameBase + '_GiveUpRate.png');
axis equal
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end
cbr.TickLabels = compose('%.1f',cbr.Ticks);

exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);

%% Plot initial food density
datafilename = filenameBase + '_FoodDensity';
psi = readFromFile([ny nx], precision, datafilename);

psimin = min(psi.*obstacle, [], 'all');
psimax = max(psi.*obstacle, [], 'all');
levels = linspace(psimin, psimax, nlevels);

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), psi((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
else
    contourf(xx, yy, psi.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2, 'LineColor', 'w');
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;
cbr = colorbar;
if psimin == psimax
    cbr.TicksMode = 'Manual';
    cbr.Ticks = psimin;
    tix = cbr.Ticks;
    cbr.TickLabels = compose('%9.2f',tix);  
end


plotfilename = strcat(outputFilenameBase + '_InitialFoodDensity.png');
axis equal
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end
cbr.TickLabels = compose('%.1f',cbr.Ticks);

exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);

%% Plot food density after stage
stageNum = 0;
datafilename = filenameBase + '_Stage' + string(stageNum) + '_FoodDensity';
psi = readFromFile([ny nx], precision, datafilename);

psimin = min(psi.*obstacle, [], 'all');
psimax = max(psi.*obstacle, [], 'all');
levels = linspace(psimin, psimax, nlevels);

homefilename = filenameBase + '_HomeBase';
homeBase = readFromFile([ny nx], precision, homefilename);

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), psi((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
else
    contourf(xx, yy, psi.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2, 'LineColor', 'w');
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;
cbr = colorbar;
if psimin == psimax
    cbr.TicksMode = 'Manual';
    cbr.Ticks = psimin;
end

plotfilename = strcat(outputFilenameBase + '_Stage' + string(stageNum) ...
             + '_FoodDensity.png');
axis equal
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end
exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);

%% Plot initial food accumulated
datafilename = filenameBase + '_FoodAccumulated';
F = readFromFile([ny nx], precision, datafilename);
Fmin = min(F.*obstacle, [], 'all');
Fmax = max(F.*obstacle, [], 'all');
levels = linspace(Fmin, Fmax, nlevels);

figure()
if example == "RealWorld"
    % Remove blank space for real world evironment
    nYShift = round(yShift*(ny-1));
    contourf(xx, yy(1:(end-nYShift)), F((nYShift+1):end,:)...
             .*obstacleMask((nYShift+1):end,:), nlevels); hold on;
    contour(xx, yy(1:(end-nYShift)), homeBase((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'w');
    contour(xx, yy(1:(end-nYShift)), obstacle((nYShift+1):end,:), [1 1], ...
            'LineWidth', 2, 'LineColor', 'k');
elseif example == "RiskReward"
    if filled_color
        contourf(xx, yy, F.*obstacle, nlevels); hold on;
        contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.2, 'LineColor', 'w');
    else
        contour(xx, yy, F.*obstacle, nlevels-2, 'LineWidth', 2); hold on;
        contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.2, 'LineColor', WongBlue);
    end
elseif example == "MultiStage"
    contourf(xx, yy, F.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.75, 'LineColor', 'k');
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2.5, 'LineColor', 'w');
    yticks([0 0.2 0.4 0.6])
else
    contourf(xx, yy, F.*obstacle, levels); hold on;
    contour(xx, yy, homeBase, [1 1], 'LineWidth', 2, 'LineColor', 'w');
end

% Figure formatting
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;
cbr = colorbar;
clim([Fmin round(Fmax)]);

plotfilename = strcat(outputFilenameBase + '_FoodAccumulated.png');
axis equal
if example == "RealWorld"
    xlim([0, 1]);
    ylim([0,0.7]);
    xticks([0 0.5 1]);
    yticks([0 0.35 0.7]);
    set(gcf,'Position', [100,100,1000,750])
else
    set(gcf,'Position', [100,100,400,300])
end
cbr.TickLabels = compose('%.1f',cbr.Ticks);

exportgraphics(gcf, plotfilename, 'Resolution', resolutionDPI);