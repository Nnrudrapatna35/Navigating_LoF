% Initialization
initialization;
pathFolder = folder + "/RR_FoodDepletion";

% Choose focal risk premium value
riskPremium = 2.0;
riskPremiumIndex = 10*riskPremium;
pathParams = "_9pt_UpperRateShift_" + string(riskPremiumIndex);

% Generate list of initial energy levels
numE0 = 20;
e0List = linspace(minE, maxE, round(numE0)+1);
e0Indices = 1:1:numE0;

% Number of paths generated
numPaths = 100;

finalUtilities = zeros(numE0, numPaths);
finalUtilitiesNoDep = zeros(numE0, numPaths);
figure;
hold on;
for k = 1:numPaths
    for j = 1:(length(e0Indices))
        e0 = e0Indices(j);
        currentParams = pathParams + "_e0_" + e0;
        pathFilenameBase = pathFolder + '/Example' + string(example) + '_'...
                           + objective + currentParams;
        
        % Read in final energy with food depletion
        energyfile = pathFilenameBase + "_Path_" + string(k) + "_Energy";
        stepsfile = pathFilenameBase + "_Path_" + string(k) + "_Steps";
    
        nSteps = readFromFile(1, "int", stepsfile);
        energy = readFromFile([nSteps, 1], precision, energyfile);
        finalUtilities(j,k) = utility(energy(nSteps));

        % Read in final energy without food depletion
        energyfile = pathFilenameBase + "_NoDepletion_Path_" + string(k) + "_Energy";
        stepsfile = pathFilenameBase + "_NoDepletion_Path_" + string(k) + "_Steps";
    
        nSteps = readFromFile(1, "int", stepsfile);
        energy = readFromFile([nSteps, 1], precision, energyfile);
        finalUtilitiesNoDep(j,k) = utility(energy(nSteps));
    end
    scatter(e0List(2:end)/maxE, finalUtilities(:, k), 50, 'filled', ...
        MarkerFaceColor=WongRedPurple, ...
        MarkerEdgeColor=WongRedPurple, MarkerFaceAlpha = 0.6, ...
        MarkerEdgeAlpha = 1.0, Marker='o', LineWidth=1);
    scatter(e0List(2:end)/maxE, finalUtilitiesNoDep(:,k), 50, 'filled', ...
        MarkerFaceColor=WongSkyBlue, ...
        MarkerEdgeColor=WongSkyBlue, MarkerFaceAlpha = 0.6, ...
        MarkerEdgeAlpha = 1.0, Marker='d', LineWidth=1);
end

meanUtility = mean(finalUtilities,2);
meanUtilityNoDep = mean(finalUtilitiesNoDep,2);
plot(e0List(2:end)/maxE, meanUtility, Color=WongRedPurple, LineWidth=2);
plot(e0List(2:end)/maxE, meanUtilityNoDep, Color=WongSkyBlue, LineWidth=2);

figs = get(gca, 'Children');
lgd = legend([figs(1), figs(2)], 'Without depletion', 'With depeletion');
lgd.Location = 'best';
xlabel('$e_0/E$', 'Interpreter','latex', 'Fontsize', labelfontsize)
ylabel('Utility', 'Interpreter','latex', 'Fontsize', labelfontsize)

ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = fontname;

if strcmp(objective, "Linear") == 1
    titles = sprintf('Risk-neutral');
elseif strcmp(objective, "Sigmoid") == 1
    titles = sprintf('Sigmoid');
elseif strcmp(objective, "SquareRoot") == 1
    titles = sprintf('Risk-averse');
end
title(titles, 'FontSize', titlefontsize, 'FontName', fontname)

set(gcf,'Position', [100,100,400*1.3,300*1.3])
saveas(gcf, outputFilenameBase + "_DepletionComparison.png");
