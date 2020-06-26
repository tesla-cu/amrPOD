clc; clear; close all;

% AMR information
dim    = 3;
finest = 1;

% Repetition information
c_l    = zeros(finest+1, 1);
for i = 0:finest
    c_l(i+1) = 2^(finest - i);
end

% =========================================================================
% Two dimensional reshaping
% =========================================================================
if dim == 2
    % Initialize grid
    grid = zeros(8,4);
    nx = size(grid, 1);
    ny = size(grid, 2);
    grid(1:2, 3:4) = 1;
    grid(5:6, 3:4) = 2;
    
    % Reshape into column vector
    for i = 1 : numel(c_l)
        nxr  = size(grid, 1);
        nyr  = size(grid, 2);
        c    = c_l(i);
        grid = permute(grid, [2, 1]);
        grid = reshape(grid, nyr, c, []);
        grid = permute(grid, [1 3 2]);
        grid = reshape(grid, [], c);
    end
    
    % Save column vector
    grid_1D = grid;
    
    % Reshape into original grid
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
    grid = reshape(grid, ny, [], nyr);
    grid = permute(grid, [1 3 2]);
    grid = reshape(grid, ny, nx);
    grid = permute(grid, [2 1]);

% =========================================================================
% Three dimensional reshaping
% =========================================================================
else  
    % Initialize grid
    grid = zeros(4,8,8);
    nx = size(grid, 1);
    ny = size(grid, 2);
    nz = size(grid, 3);
    grid(1:2, 1:2, 1:2) = 2;
    grid(3:4, 1:2, 1:2) = 1;
    nspat = numel(grid);
    
    % Reshape into column vector
    for i = 1 : numel(c_l)
        c       = c_l(i);
        nxr     = size(grid, 1);
        nyr     = size(grid, 2);
        nzr     = size(grid, 3);
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
    
    % Save column vector
    grid_1D = grid;
    
    % Reshape into original grid
    for i = numel(c_l)-1 : -1 : 1 
        c       = c_l(i);
        nxr     = size(grid, 1);
        nyr     = size(grid, 2);
        nzr     = size(grid, 3);
        nrshp_1 = nspat / (c*nyr*nzr);
        nrshp_2 = nspat / (c*c*nzr);
        nrshp_3 = nspat / (c*c);
        
        grid = reshape(grid, c, nrshp_1, nyr, nzr);
        grid = permute(grid, [1 3 2 4]);
        grid = reshape(grid, c, c, nrshp_2, nzr);
        grid = permute(grid, [1 2 4 3]);
        grid = reshape(grid, c, c, nrshp_3);
        grid = permute(grid, [3 2 1]);
    end
    
    nxr     = size(grid, 1);
    nyr     = size(grid, 2);
    nzr     = size(grid, 3);
    nrshp_1 = nspat/(nz*nyr*nzr);
    nrshp_2 = nspat/(nz*ny*nzr);
    
    grid = reshape(grid, nz, nrshp_1, nyr, nzr);
    grid = permute(grid, [1 3 2 4]);
    grid = reshape(grid, nz, ny, nrshp_2, nzr);
    grid = permute(grid, [1 2 4 3]);
    grid = reshape(grid, nz, ny, nx);
    grid = permute(grid, [3 2 1]);
    
end