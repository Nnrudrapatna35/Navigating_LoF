clear all; close all; clc;
%% Generate Obstacle Data
% Note - trying to recreate what was done during the REU, this doesn't work
% yet
filename = "Food.png"; % Has a clear black outline
colorData = imread(filename);
imshow(colorData)
colorData = double(colorData);
colorData = flip(colorData);

[nx, ny, nc] = size(colorData);
N = max([nx, ny]);
obstacleXIndex = [];
obstacleYIndex = [];

% R = colorData(:,:,1); 
% figure
% contourf(R)
% colorbar
% G = colorData(:,:,2);
% figure
% contourf(G)
% colorbar
% B = colorData(:,:,3);
% figure
% contourf(B)
% colorbar

obstacle = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        rgbData = colorData(i,j,:);
        rgbNorm = norm(reshape(rgbData, [],3));

        % Colors near black will have a small norm, colors near white will
        % have a large norm, and colors near grey will have similar r, g,
        % and b values. Thresholds based on trial and error.
        if rgbNorm < 140
            obstacle(i,j) = 1;
        elseif rgbNorm > 410
            obstacle(i,j) = 1;
        elseif (abs(rgbData(1)-rgbData(2)) + abs(rgbData(2)-rgbData(3)) ...
                + abs(rgbData(3)-rgbData(1))) < 70
            obstacle(i,j) = 1;
        end
    end
end

obstacle = 1 - obstacle;
paddingX = zeros(N-nx, ny);
paddingY = zeros(nx, N-ny);

obstacle = [obstacle paddingY];
obstacle = [paddingX; obstacle];

figure
contourf(obstacle)
writematrix(obstacle, "Obstacle.csv");
%% Generate Psi data (given obstacle data)
% Read in image data
filename = "Food.png";
F = imread(filename);

% Read in obstacle data
obstacle = readmatrix('Obstacle.csv');
isObstacle = (obstacle==0);

% Get RGB data
F = double(F);
R = F(:,:,1); 
G = F(:,:,2);
B = F(:,:,3);

% Compute Psi as red (high) - green (low) and rescale
psi = R - G;
psiMax = max(psi, [], 'all'); 
psiMin = min(psi, [], 'all');
psi = (1/(psiMax-psiMin)).*(psi-psiMin);

% Interpolate to desired resolution and make the domain square
% [nx, ny] = size(psi);
% [xx, yy] = ndgrid(1:nx,1:ny);
% psiInterpolator = griddedInterpolant(xx,yy,psi);
% 
% xq = 0:0.5:nx; % We choose to upscale the resolution
% yq = 0:0.5:ny;
% [xxq, yyq] = ndgrid(xq, yq);
% psiRefined = psiInterpolator(xxq, yyq);

% Fix orientation, shape, and zero out at obstacle
% psi((2*(nx+1)):(2*ny+1),:)=0;
psi=flip(psi);
psi = [psi paddingY];
psi = [paddingX; psi];

psi(isObstacle)=0;

figure
contourf(psi)
writematrix(psi, "Psi.csv");

%% Generate MuS data (given obstacle data)
% Read in image data
filename = "Fear.png";
F = imread(filename);
F = imgaussfilt(F,2);

% Read in obstacle data
obstacle = readmatrix('Obstacle.csv');
isObstacle = (obstacle==0);

% Get RGB data
F = double(F);
R = F(:,:,1); 
G = F(:,:,2);
B = F(:,:,3);

% Compute Psi as red (high) - blue (low) and rescale
mus = R - B;
musMax = max(mus, [], 'all'); 
musMin = min(mus, [], 'all');
mus = (1/(musMax-musMin)).*(mus-musMin);

% Interpolate to desired resolution and make the domain square
% [nx, ny] = size(mus);
% [xx, yy] = ndgrid(1:nx,1:ny);
% musInterpolator = griddedInterpolant(xx,yy,mus);
% 
% xq = 0:0.5:nx; % We choose to upscale the resolution
% yq = 0:0.5:ny;
% [xxq, yyq] = ndgrid(xq, yq);
% musRefined = musInterpolator(xxq, yyq);

% Fix orientation, shape, and zero out at obstacle
% musRefined((2*(nx+1)):(2*ny+1),:)=0;
mus=flip(mus);
mus = [mus paddingY];
mus = [paddingX; mus];
mus(isObstacle)=0;

figure
surf(mus);
writematrix(mus, "MuS_smooth.csv");

%% Compute MuG Data
% mus = readmatrix('MuS.csv');
obstacle = readmatrix('Obstacle.csv');
isObstacle = (obstacle==0);

mug = 1 - mus;
mug(isObstacle) = 0;

figure
contourf(mug);
writematrix(mug, "MuG.csv");

%% Sleeping sites
% Note - trying to recreate what was done during the REU, this doesn't work
% yet
filename = "SleepingSites.png"; % Has a clear black outline
colorData = imread(filename);
colorData = double(colorData);
colorData = flip(colorData);

obstacle = readmatrix('Obstacle.csv');
isObstacle = (obstacle==0);

[nx, ny, nc] = size(colorData);
N = max([nx, ny]);

R = colorData(:,:,1);
G = colorData(:,:,2);
B = colorData(:,:,3);

sites = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        rgbData = colorData(i,j,:);
        rgbNorm = norm(reshape(rgbData, [],3));

        % Want to isolate areas that are blue
        if B(i,j) > 130 && R(i,j) < 170
            sites(i,j) = 1;
        end
    end
end

% sites = 1 - sites;
sites = [sites paddingY];
sites = [paddingX; sites];
sites(isObstacle)=0;

figure
contourf(sites)
writematrix(sites, "SleepingSites.csv");

%% Scratchwork
% Get name of image
filename="Food.jpg";

% Read it into a RGB matrix
F = imread(filename);
figure;
image(F)

H = rgb2hsv(F);
H = flip(H);
F = double(F);
% F = flip(F);
sz = size(F);
X = 1:sz(1);
Y = 1:sz(2);


FRed = F(:,:,1);
FGreen = F(:,:,2);
FBlue = F(:,:,3);

figure
subplot(1,3,1);
contourf(FRed);
title('Red')
colorbar;
subplot(1,3,2);
contourf(FGreen);
title('Green')
colorbar;
subplot(1,3,3);
contourf(FBlue);
colorbar;
title('Blue')

HHue = H(:,:,1);
HSat = H(:,:,2);
HVal = H(:,:,3);

figure
subplot(1,3,1);
contourf(HHue);
colorbar;
title('Hue')
subplot(1,3,2);
contourf(HSat);
colorbar;
title('Saturation')
% contour(HSat, [-0.1 0.1]) % Good option for finding obstacle
subplot(1,3,3);
contourf(HVal);
colorbar;
title('Value')
