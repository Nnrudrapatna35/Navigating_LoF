% Plots the optimal patch behavior as a function of initial energy and risk
% premium for the Risk Reward example. Choose "paper" for a larger range of
% parameter values and "poster" for a smaller range of parameters

initialization;
pathFolder = folder + "/RR_PatchSelection";

% Choose starting, ending, and increment risk premium values
minRiskPremium = 0.1;
maxRiskPremium = 3.0;
riskPremiumIncrement = 0.1;

% Convert risk premium to integer labels
riskPremiumIndices = 10*(minRiskPremium:riskPremiumIncrement:maxRiskPremium);

% Generate list of initial energy levels
numE0 = 20;
e0List = linspace(minE, maxE, round(numE0)+1);
e0Indices = 0:10*riskPremiumIncrement:numE0;

% Choose energy levels to display
minDisplayEnergy = 3.0;
e0Indices = e0Indices(e0List>=minDisplayEnergy);
e0List = e0List(e0List>=minDisplayEnergy)/maxE;

% Store patch choice
patchChoice = zeros(length(riskPremiumIndices), length(e0List));
for i= 1:length(riskPremiumIndices)
    riskPremium = riskPremiumIndices(i);
    for j = 1:(length(e0List))
        e0 = e0Indices(j);

        currentParams = params + "_UpperRateShift_" ...
                      + string(riskPremium) + "_e0_" + e0 ...
                      + "_NoDepletion";
        pathFilenameBase = pathFolder + "/Example" + string(example)...
                            + "_" + objective + currentParams;
        
        typefile = pathFilenameBase + "_PathType";
        stepsfile = pathFilenameBase + "_Steps";

        patchChoice(i,j) = readFromFile(1, "int", typefile);
        nSteps = readFromFile(1, "int", stepsfile);
    end
end

% Plot results
figure();
[xx, yy] = ndgrid(0.1*riskPremiumIndices, e0List);
xAxisValues = reshape(xx.', 1, []); 
yAxisValues = reshape(yy.', 1, []); 
pathChoice = reshape(patchChoice.', 1, []);

% Plot each patch selection behavior
symbolSize = 1.2;
hold on;

% Upper patch
xUpper = xAxisValues(pathChoice == 1);
yUpper = yAxisValues(pathChoice == 1);
scatter(xUpper, yUpper, [], "^", 'filled', 'MarkerFaceColor',...
        (WongVermillion + [0.8 0.1 0.3])/2, 'MarkerEdgeColor', ...
        (WongVermillion + [0.8 0.1 0.3])/2, 'LineWidth', symbolSize); 

% Lower patch
xLower = xAxisValues(pathChoice == 2);
yLower = yAxisValues(pathChoice == 2);
scatter(xLower, yLower, [], "v",'filled', 'MarkerFaceColor', WongOrange,...
        'MarkerEdgeColor', WongOrange, 'LineWidth',symbolSize); 

% Stay at home
xHome = xAxisValues(pathChoice == 0);
yHome = yAxisValues(pathChoice == 0);
scatter(xHome, yHome, [], "o",'filled', 'MarkerFaceColor', WongBlueGreen,...
        'MarkerEdgeColor', WongBlueGreen, 'LineWidth',symbolSize); 

% Visit both patches
xBoth = xAxisValues(pathChoice == 3);
yBoth = yAxisValues(pathChoice == 3);
scatter(xBoth, yBoth, [], "*",'filled', 'MarkerFaceColor', WongYellow,...
        'MarkerEdgeColor', WongYellow, 'LineWidth',symbolSize); 
hold off;

% Figure formatting
ylim([minDisplayEnergy/maxE - 0.05, maxE + 0.05]) 
xlim([minRiskPremium - riskPremiumIncrement, maxRiskPremium ...
      + riskPremiumIncrement])
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = 'LM Roman 10';

xlabel('Risk Premium in Upper Half');
h = ylabel('$e_{0}/E$', 'Interpreter','latex');
h.FontSize = labelfontsize;

set(gcf,'Position', [100,100,500*1.2,300*1.2])
plotfilename = strcat(outputFilenameBase + '_EnergyThresholdPlot');
saveas(gcf, plotfilename + ".png");