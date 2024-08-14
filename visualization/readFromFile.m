function grid = readFromFile(dimensions, precision, filename)
% Order of dimensions: [ne, ny, nx, nt] if 4, [ne, ny, nx] if 3, [ny, nx]
% if 2
numDim = size(dimensions);

file = fopen(filename);
grid = fread(file, prod(dimensions), precision); 
fclose(file);

if numDim(2) >= 2
    grid = reshape(grid, dimensions);
end

if numDim(2) == 4
    grid = permute(grid, [2 3 1 4]);
elseif numDim(2) == 3
    grid = permute(grid,[2 3 1]);
end
end