%% Graph Utility functions
initialization

eThreshold = 0.3;
maxEnergy = 1;
theta = 0.1;
ee = linspace(0, maxEnergy, 101);

figure;
utility = @(e) e;
plot(ee, utility(ee), "Color", WongYellow, "LineWidth", 2.5); hold on;

utility = @(e) sqrt(e);
plot(ee, utility(ee), "Color", WongBlueGreen, "LineWidth", 2.5, "LineStyle","--");

utility = @(e) (erf((e - eThreshold)/(sqrt(2)*theta)) - erf((-eThreshold)/(sqrt(2)*theta)))...
          / (erf((maxEnergy - eThreshold)/(sqrt(2)*theta)) - erf((-eThreshold)/(sqrt(2)*theta)));
plot(ee, utility(ee), "Color", WongSkyBlue, "LineWidth", 2.5, "LineStyle","-.");

legend("Risk-neutral", "Risk-averse", "Sigmoid", "Location","southeast")
xlabel("$e$", "Interpreter","latex")
ylabel("$U(e)$", "Interpreter","latex")
xlim([0,1])
ax = gca;
ax.FontSize = labelfontsize;
ax.FontName = figurefont;

% axis equal
set(gcf,"Position",[100 100 400*(1) 300*(1)])
exportgraphics(gcf, "../plots/UtilityFunctions.eps")%, 'Resolution', 900);

