
clc;clear;
dim    = 3;
finest = 1;
c_l    = zeros(finest+1, 1);
d_l    = zeros(finest+1, 1);
for i = 0:finest
    c_l(i+1) = 2^(finest - i);
    d_l(i+1) = (2^dim)^(finest-i);
end

if dim == 2
    grid = zeros(8,4);
    nx = size(grid, 1);
    ny = size(grid, 2);
    grid(1:2, 3:4) = 1;
    grid(5:6, 3:4) = 2;
    for i = 1 : numel(c_l)
        nxr  = size(grid, 1);
        nyr  = size(grid, 2);
        c    = c_l(i);
        grid = permute(grid, [2, 1]);
        grid = reshape(grid, nyr, c, []);
        grid = permute(grid, [1 3 2]);
        grid = reshape(grid, [], c);
    end
    
    grid_1D = grid;
    
    for i = numel(c_l)-1 : -1 : 1 
        nxr  = size(grid, 1);
        nyr  = size(grid, 2);
        c    = c_l(i);
        grid = reshape(grid, c, [], nyr);
        grid = permute(grid, [1 3 2]);
        grid = reshape(grid, c, []);
        grid = permute(grid, [2 1]);
    end
    nxr  = size(grid, 1);
    nyr  = size(grid, 2);
    c = c(1);
    grid = reshape(grid, ny, [], nyr);
    grid = permute(grid, [1 3 2]);
    grid = reshape(grid, ny, nx);
    grid = permute(grid, [2 1]);
    
else
    grid = zeros(4,8,8);
    nx = size(grid, 1);
    ny = size(grid, 2);
    nz = size(grid, 3);
    grid(1:2, 1:2, 1:2) = 2;
    grid(3:4, 1:2, 1:2) = 1;
    nspat = numel(grid);
    
    for i = 1 : numel(c_l)
        c = c_l(i);
        nsize   = size(grid);
        nxr     = nsize(1);
        nyr     = nsize(2);
        nzr     = nsize(3);
        nrshp_1 = nspat / (nzr*nyr*c);
        nrshp_2 = nspat / (nzr*c^2);
        nrshp_3 = nspat / c^2;
        
        grid = permute(grid, [3, 2, 1]);
        grid = reshape(grid, nzr, nyr, c, nrshp_1);
        grid = permute(grid, [1, 2, 4, 3]);
        grid = reshape(grid, nzr, c, nrshp_2, c);
        grid = permute(grid, [1, 3, 2, 4]);
        grid = reshape(grid, nrshp_3, c, c);
    end
    
    grid_1D = grid;
    
    for i = numel(c_l)-1 : -1 : 1 
        nxr  = size(grid, 1);
        nyr  = size(grid, 2);
        nzr  = size(grid, 3);
        c    = c_l(i);
        grid = reshape(grid, c, [], nyr, nzr);
        grid = permute(grid, [1 3 2 4]);
        grid = reshape(grid, c, c, [], nzr);
        grid = permute(grid, [1 2 4 3]);
        grid = reshape(grid, c, c, []);
        grid = permute(grid, [3 2 1]);
    end
    nxr  = size(grid, 1);
    nyr  = size(grid, 2);
    nzr  = size(grid, 3);
    c = c(1);
    grid = reshape(grid, nz, [], nyr, nzr);
    grid = permute(grid, [1 3 2 4]);
    grid = reshape(grid, nz, ny, [], nzr);
    grid = permute(grid, [1 2 4 3]);
    grid = reshape(grid, nz, ny, nx);
    grid = permute(grid, [3 2 1]);
    
end