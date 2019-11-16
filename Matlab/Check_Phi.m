clc; clear; close all;

% Grid information
nx    = 512;
ny    = 512;
nz    = 1;
nspat = nx*ny*nz;
% nt    = 10;
nt    = 101;

% Direction where /code/ lives
basedir = '/Users/mikemeehan/Research/Papers/2019_POD_AMR/';

% Directory where AMR data is stored
amr_datadir = ['/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/data/' ...
    'slice/x0.000-0.000_y-1.000-1.000_z0.000-2.000_t40.0000-42.0000/'];

% Directory where we was to store data on speed up
datadir = [basedir 'data/'];

% Directory that describes the study we are looking at
studydir = [datadir 'check_MAT/'];

% File information
file_fmt = '%s%05i.bin';
var      = 'z_velocity';

% Load data
X = zeros(nspat, nt);
for i = 0 : nt-1
    file = sprintf([amr_datadir '%s%05i.bin'],var, i);
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

% Load python data
R_py   = load([studydir 'R_py.txt']);
Phi_py = load([studydir 'Phi_py.txt']);
A_py   = load([studydir 'A_py.txt']);

% Compare differences between Matlab and python (note, need two sets of abs
% for Phi and A because these can be identically opposite, which is still
% correct)
R_diff   = max(max(abs(R   - R_py)));
Phi_diff = max(max(abs(Phi) - abs(Phi_py)));
A_diff   = max(max(abs(A)   - abs(A_py)));
fprintf(['The maximum difference between the matrices computed with python ' ...
    'and Matlab is:' ...
    '\n\t R:   %0.8e' ...
    '\n\t Phi: %0.8e' ...
    '\n\t A:   %0.8e\n'], R_diff, Phi_diff, A_diff);



