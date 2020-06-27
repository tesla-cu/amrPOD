clc; clear; close all;

% Grid information
nx    = 128;
ny    = 128;
nz    = 128;
nspat = nx*ny*nz;
nt    = 5;

% Direction where /code/ lives
basedir = '../../../';

% Directory where AMR data is stored
amr_datadir = ['/Users/mikemeehan/Research/Papers/2019_POD_AMR/data/'...
    'AMR_gen/f1_3D/'];

% Directory where POD data computed by fortran code is stored
PODdir = [amr_datadir 'POD_data_nt5/'];

% Directory where we was to store data on speed up
datadir = [basedir 'data/'];

% Directory that describes the study we are looking at
studydir = [datadir 'fortran_vs_matlab/'];

% File information
file_fmt = '%s%05i.bin';
var      = 'z_velocity';

% Load data
X = zeros(nspat, nt);
for i = 0 : nt-1
    file = sprintf([amr_datadir file_fmt],var, i);
    fid = fopen(file);
    var_data = fread(fid, nspat, 'float64');
    X(:,i+1) = var_data;
    fclose(fid);
end

% Compute POD
R            = X'*X;
[Psi,Lambda] = eig(R);
Phi          = X*Psi/sqrt(Lambda);
A            = X'*Phi;

% Load Phi fortran data
Phi_fort = zeros(nspat, nt);
for i = 0 : nt-1
    file = sprintf([PODdir file_fmt], ['POD_' var], i);
    fid = fopen(file);
    var_data = fread(fid, nspat, 'float64');
    Phi_fort(:,i+1) = var_data;
    fclose(fid);
end

% Load A fortran data
file = [PODdir 'temporal_coefficients.bin'];
fid = fopen(file);
A_fort = fread(fid, nt*nt, 'float64');
A_fort = reshape(A_fort, [nt,nt]);
fclose(fid);

% Compare differences between Matlab and Fortran (note, need two sets of abs
% for Phi and A because these can be identically opposite, which is still
% correct)
Phi_diff = max(max(abs(Phi) - abs(Phi_fort)));
A_diff   = max(max(abs(A)   - abs(A_fort)));
fprintf(['The maximum difference between the matrices computed with ' ...
    'fortran and Matlab is:' ...
    '\n\t Phi: %0.8e' ...
    '\n\t A:   %0.8e\n'], Phi_diff, A_diff);



