clearvars; close all; clc;

% Set names and parmeters for the example
example = "MultiStage"; % specify: RiskReward, RiskRewardPT, RealWorld, MultiStage
stageNum = 0; %if single stage, specify stage number
modeNum = 1;
objective = "Linear";
if strcmp(objective, "Linear") == 1
    objectiveName = 'Risk-neutral';
elseif strcmp(objective, "Sigmoid") == 1
    objectiveName = 'Sigmoid';
elseif strcmp(objective, "SquareRoot") == 1
    objectiveName = 'Risk-averse';
elseif strcmp(objective, "Utility") == 1
    objectiveName = 'Utility';
end

folder = "../output";
params = "_9pt";
yShift = 0;
if example == "RiskRewardPT"
    params = "_9pt_RateScaled_20";
    example = "RiskReward";
elseif example == "RealWorld"
    % For real world example
    yShift = 0.3;
elseif example == "MultiStage"
    params = "_9pt_N_5_mus_3_muk_3_realization_0";
end

filenameBase = folder + '/Example' + string(example) + '_' + objective + params;
precision = "double";
valuePrecision = "float";

% Figure settings
set(0, 'DefaultFigureRenderer', 'painters')
titlefontsize = 20;
labelfontsize = 16;
largelabelfontsize = 32;
figurefont = 'LM Roman 10';
resolutionDPI = 600;

% Create output directory if needed
outputParams = params;
if modeNum == 2
    outputParams = params + "_Mode2";
end
outputFolder = "Example" + string(example);
if ~exist("../plots/" + outputFolder, 'dir')
    mkdir("../plots/" + outputFolder)
end
outputFilenameBase = "../plots/" + outputFolder + '/Example' ...
    + string(example) + '_' + objective + outputParams;

% Important input filenames
paramsfile = filenameBase + '_Stage' + string(stageNum) + '_Mode1Parameters';
initialfoodfile = filenameBase + '_Stage' + string(stageNum) + '_FoodDensity';
runningcostfile = filenameBase + "_Stage" + string(stageNum) + "_Mode1Cost";
speedfile = filenameBase + "_Stage" + string(stageNum) + "_Mode1Speed";
timefactorfile = filenameBase + "_Stage" + string(stageNum) + "_TimeFactor";

% Read in parameters
file = fopen(paramsfile);
p = fread(file, precision);
fclose(file);

nx = p(1);
ny = p(2);
ne = p(3);
nt = p(4);
minXY = p(5);
maxXY = p(6);
dt = p(7);
dx = p(8); 
minE = p(9);
maxE = p(10);
de = (maxE - minE)/(ne-1);

% Define color palette
WongBlack = [0/255   0/255   0/255]; % Black
WongOrange = [230/255 159/255   0/255]; % Orange
WongSkyBlue = [86/255 180/255 233/255]; % Sky blue
WongBlueGreen = [0/255 158/255 115/255]; % Bluish green
WongYellow = [240/255 228/255  66/255]; % Yellow
WongBlue = [0/255 114/255 178/255]; % Blue
WongVermillion = [213/255  94/255   0/255]; % Vermillion
WongRedPurple = [204/255 121/255 167/255]; % Reddish purple

% Utility Functions
eThreshold = 0.3*maxE;
theta = 0.1*maxE;
if strcmp(example, "9")
    theta = 0.1*maxE;
    eThreshold = 0.5*maxE;
end

utilityM = @(e) (erf((e - eThreshold)/(sqrt(2)*theta)) - erf((-eThreshold)/(sqrt(2)*theta)))...
          / (erf((maxE - eThreshold)/(sqrt(2)*theta)) - erf((-eThreshold)/(sqrt(2)*theta)));
utilityR = @(e) sqrt(e/maxE);
utilityL = @(e) e/maxE;
utilityE= @(e) e;

if strcmp(objective, "SquareRoot")
    utility = @(e) utilityR(e);
elseif strcmp(objective, "Sigmoid")
    utility = @(e) utilityM(e);
elseif strcmp(objective, "Linear")
    utility = @(e) utilityL(e);
elseif strcmp(objective, "Energy")
    utility = @(e) utilityE(e);
elseif strcmp(objective, "Utility")
    eThreshold = 0.5*maxE;
    utility = @(e) (erf((e - eThreshold)/(sqrt(2)*theta)) - erf((-eThreshold)/(sqrt(2)*theta)))...
          / (erf((maxE - eThreshold)/(sqrt(2)*theta)) - erf((-eThreshold)/(sqrt(2)*theta)));
end

nYShift = round(yShift/dx);