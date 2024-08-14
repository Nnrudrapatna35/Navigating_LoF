function value = interpolate2D(x,y,dx,dy,grid)
% 2D interpolation using the values x and y, assuming points in grid are dx
% apart in the first dimension and dy apart in the second dimension. 
% Assumes x and y are measured from 0.

numPoints = size(grid);
iLow = floor(x/dx) + 1;
jLow = floor(y/dy) + 1;

if iLow == numPoints(1)
    iLow = iLow - 1;
end

if jLow == numPoints(2)
    jLow = jLow - 1;
end

Cx = (x - dx*(iLow-1))/dx;
Cy = (y - dy*(jLow-1))/dy;

value = (1 - Cx)*((1 - Cy)*grid(iLow, jLow) + Cy*grid(iLow, jLow + 1))...
        + Cx*((1 - Cy)*grid(iLow + 1, jLow) + Cy*grid(iLow+1, jLow+1));
end