clc; clear; close all;

% =========================================================================
% User defined inputs
% =========================================================================

% Set finest level on grid
finest = 2;

% Two dimensional grid (comment out if using 3D grid)
grid = zeros(8,8); 
grid(1:2, 3:4) = 1.1 + rand(2,2)/1000;
grid(5:6, 3:4) = 1.2 + rand(2,2)/1000;
grid(3:4, 1:2) = 1.3 + rand(2,2)/1000;

% Three dimensional grid (comment out if using 2D grid)
% grid = zeros(8,8,8);
% grid(1:2, 1:2, 1:2) = 1.1 + rand(2,2,2)/100;
% grid(3:4, 1:2, 1:2) = 1.2 + rand(2,2,2)/100;
% grid(1:2, 1:2, 3:4) = 1.3 + rand(2,2,2)/100;
% grid(5:8, 5:8, 5:8) = 0.1 + rand(4,4,4)/100;

% =========================================================================
% Define useful information
% =========================================================================

% Dimension of data
dim = length(size(grid));

% Store original version of the grid
grid_orig = grid;

% Compute number of cells in each directon
nx    = size(grid, 1);
ny    = size(grid, 2);
nspat = nx*ny;
if dim == 3; nz = size(grid, 3); nspat = nspat*nz; end

% Repetition information
c_l = zeros(finest+1, 1);
for i = 0:finest
    c_l(i+1) = 2^(finest - i);
end

% =========================================================================
% Two dimensional reshaping
% =========================================================================
if dim == 2
    
    % Reshape into column vector ------------------------------------------
    for i = 1:length(c_l)
        c     = c_l(i);          % c value for iteration number
        nx1   = size(grid, 1);   % nx of initial grid
        ny1   = size(grid, 2);   % ny of initial grid
        nx2   = nspat/c;         % nx of new grid
        ny2   = c;               % ny of new grid
        grid2 = zeros(nx2, ny2); % initialize new grid
        
        % Begin assignment of new grid
        for j2=1:c
            i2 = 1;
            for j1=j2:c:ny1
                for ii=1:c:nx1
                    for i1=ii:ii+c-1
                        grid2(i2,j2) = grid(i1,j1);
                        i2 = i2 + 1;
                    end
                    i2 = i2 + ny1 - c;
                end
                i2 = i2 - nx2 + c;
            end
        end
        
        % Reassign grid for next iteration
        grid = grid2;
    end
    
    % Save column vector --------------------------------------------------
    grid_1D = grid;

    % Reshape into original grid ------------------------------------------
    for i = length(c_l)-1:-1:0
        nx1  = size(grid, 1);    % nx of initial grid
        ny1  = size(grid, 2);    % ny of initial grid
        % Check if we are on the final iteration
        if i > 0 
            c    = c_l(i);       % c value for iteration number
            nx2 = nspat/c;       % nx of new grid
            ny2 = c;             % ny of new grid
        else
            nx2 = nx;            % nx of original grid
            ny2 = ny;            % ny of original grid
        end
        grid2 = zeros(nx2, ny2); % initialize new grid
        
        % Begin assignment of new grid
        for j1=1:ny1
            i1 = 1;
            for j2=j1:ny1:ny2
                for ii=1:ny1:nx2
                    for i2=ii:ii+ny1-1
                        grid2(i2,j2) = grid(i1,j1);
                        i1 = i1 + 1;
                    end
                    i1 = i1 + ny2 - ny1;
                end
                i1 = i1 - nx1 + ny1;
            end
        end
        
        % Reassign grid for next iteration
        grid = grid2;
    end
end
    
% =========================================================================
% Three dimensional reshaping
% =========================================================================
if dim == 3  
    
    % Reshape into column vector ------------------------------------------
    for i = 1:length(c_l)
        c     = c_l(i);               % c value for iteration number     
        nx1   = size(grid, 1);        % nx of initial grid 
        ny1   = size(grid, 2);        % ny of initial grid
        nz1   = size(grid, 3);        % nz of initial grid
        nx2   = nspat/(c*c);          % nx of new grid
        ny2   = c;                    % ny of new grid
        nz2   = c;                    % nz of new grid
        grid2 = zeros(nx2, ny2, nz2); % initialize new grid
        
        % Begin assignment of new grid
        for k2=1:c
            for j2=1:c
                i2 = 1;
                for k1=k2:c:nz1
                    for j1=j2:c:ny1
                        for ii=1:c:nx1
                            for i1=ii:ii+c-1
                                grid2(i2,j2,k2) = grid(i1,j1,k1);
                                i2 = i2 + 1;
                            end
                            i2 = i2 + ny1*nz1/c - c;
                        end
                        i2 = i2 - nx2 + c;
                    end
                end
            end
        end
        
        % Reassign grid for next iteration
        grid = grid2;
    end
    
    % Save column vector --------------------------------------------------
    grid_1D = grid;
    
    % Reshape into original grid ------------------------------------------
    for i = length(c_l)-1:-1:0
        nx1  = size(grid, 1);         % nx of initial grid
        ny1  = size(grid, 2);         % ny of initial grid
        nz1  = size(grid, 2);         % nz of initial grid
        % Check if we are on the final iteration
        if i > 0
            c    = c_l(i);            % c value for iteration number
            nx2 = nspat/(c*c);        % nx of new grid
            ny2 = c;                  % ny of new grid
            nz2 = c;                  % nz of new grid
        else
            nx2 = nx;                 % nx of original grid
            ny2 = ny;                 % ny of original grid
            nz2 = nz;                 % nz of original grid
        end
        grid2 = zeros(nx2, nz2, nz2); % initialize new grid
        
        % Begin assignment of new grid
        for k1=1:nz1
            for j1=1:ny1
                i1 = 1;
                for k2=k1:nz1:nz2
                    for j2=j1:ny1:ny2
                        for ii=1:ny1:nx2
                            for i2=ii:ii+ny1-1
                                grid2(i2,j2,k2) = grid(i1,j1,k1);
                                i1 = i1 + 1;
                            end
                            i1 = i1 + ny2*nz2/ny1 - ny1;
                        end
                        i1 = i1 - nx1 + ny1;
                    end
                end
            end
        end
        
        % Reassign grid for next iteration
        grid = grid2;
    end
end


    