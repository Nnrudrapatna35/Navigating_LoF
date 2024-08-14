% Plots the optimal trajectory as a function of initial energy and for a 
% fixed risk premium for the Risk-Reward example. Can choose to show 3 or 6
% possible initial energy levels. Can show with depeletion, without
% depletion, or comparing the two.

% Initialization
initialization
pathFolder = folder + "/RR_PathDeformation";

% Compare with food depletion case?
compareDepletion = true;

% Choose focal risk premium value
riskPremium = 2.0;
riskPremiumIndex = 10*riskPremium;

% Update params
pathParams = "_9pt_UpperRateShift_" + string(riskPremiumIndex);

% Choose e0 values
e0List = [0.3 0.55 0.8];
e0Indices = 20*e0List;

% Read background
twoDimensions = [nx, ny];
mus = readFromFile(twoDimensions, precision, filenameBase + "_PredatorDensity");
F = readFromFile(twoDimensions, precision, filenameBase + "_FoodAccumulated");
home = readFromFile(twoDimensions, precision, filenameBase + "_HomeBase");

musMax = max(mus, [], 'all');
musMin = min(mus(home==1), [], 'all');

xx = linspace(minXY,maxXY,nx);
yy = linspace(minXY,maxXY,ny);

t = tiledlayout(1,3,'TileSpacing','loose','Padding','compact');
for i=1:3
    nexttile
    e0 = e0List(i);
    e0Index = e0Indices(i);
    
    % Read in path
    currentParams = pathParams + "_e0_" + e0Index;
    pathFilenameBase = pathFolder + '/Example' + string(example) ...
                     + '_' + objective + currentParams + "_NoDepletion";
    
    nSteps = readFromFile(1, "int", pathFilenameBase + "_Steps");
    path = readFromFile([2,nSteps], precision, pathFilenameBase + "_Path");
    pathX = path(1,:);
    pathY = path(2,:);
    
    % Separate "outbound" and "inbound" sections of the path
    breakPoint = find_breakpoint(pathY);

    colorout = 'w';
    colorin = 'w';
    grey = [0.5 0.5 0.5];
    symbolout = "-";
    symbolin = ":";

    % Graph results
    nlevels = 7;
    contourf(xx, yy, mus, nlevels, 'LineWidth', 0.1); hold on;
    contour(xx, yy, F, [1 1], 'LineWidth', 3, 'LineColor', 'g');
    contour(yy, xx, home, [1 1], 'LineWidth', 2, 'LineColor', 'w');

    % Plot trajectory
    plot(pathX(1:breakPoint), pathY(1:breakPoint), 'Color', colorout, ...
         'LineWidth', 3.5, 'LineStyle', symbolout);
    plot(pathX(breakPoint:end), pathY(breakPoint:end), 'Color', colorin, ...
         'LineWidth', 3.5, 'LineStyle', symbolin);
    scatter(pathX(1), pathY(1), 'o', 'filled', 'LineWidth', 1.5, ...
            'MarkerEdgeColor', grey, 'MarkerFaceColor', grey)

    % Compare with food depletion version
    if (compareDepletion)
        pathFilenameBase = pathFolder + '/Example' ...
                         + string(example) + '_' + objective ...
                         + currentParams;
        
        nSteps = readFromFile(1, "int", pathFilenameBase + "_Steps");
        path = readFromFile([2,nSteps], precision, pathFilenameBase + "_Path");
        pathX = path(1,:);
        pathY = path(2,:);
    
        breakPoint = find_breakpoint(pathY);

        plot(pathX(breakPoint:end), pathY(breakPoint:end), 'Color', 'k', ...
             'LineWidth', 3.5, 'LineStyle', "-.");
    end
    
    % Figure formatting
    title(sprintf('$e_0/E = $ %1.2f', e0), 'Interpreter','latex', ...
          'Fontsize', largelabelfontsize*0.9)
    ax = gca;
    ax.FontSize = labelfontsize*0.8;
    ax.FontName = 'LM Roman 10';
    axis equal;
end

set(gcf,'Position', [100,100,1000*(1),500*0.7])

if compareDepletion
    saveas(gcf, outputFilenameBase + "_DepletionComparison_PathDeformation.png");
else
    saveas(gcf, outputFilenameBase + "_PathDeformation_.png");
end