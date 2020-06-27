clearvars, clc
% datadir = '/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/data/volume/x-0.500-0.500_y-0.250-0.250_z0.000-1.000_t40.0000-40.0200_f2/';
% datadir = '/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/data/volume/x-0.500-0.500_y-0.250-0.250_z0.000-1.000_t40.0000-40.0200_f3/';
% datadir = '/Users/mikemeehan/Research/InternalResearchPapers/AMR_POD/data/volume/x-1.000-1.000_y-1.000-1.000_z0.000-2.000_t40.0000-40.0200_f3/';
datadir = '/Users/mikemeehan/Research/Papers/2019_POD_AMR/data/AMR_sim/';
datafile = fullfile(datadir, 'grid_level00001.bin');

% nx = 256;
% ny = 128;
% nz = 256;

% nx = 128;
% ny = 64;
% nz = 128;

nx = 512;
ny = 512;
nz = 512;

fid = fopen(datafile);
data = fread(fid, nx*ny*nz, 'float64');
data = reshape(data, [nz, ny, nx]);

data = permute(data, [3,2,1]);

figure
imagesc(squeeze(data(nx/2,:,:)))
axis equal xy tight
xlabel('z')
ylabel('y')
colorbar
% 
figure
imagesc(squeeze(data(:,ny/2,:)))
axis equal xy tight
xlabel('z')
ylabel('x')
colorbar

figure
imagesc(squeeze(data(:,:,nz/2 - 7)))
axis equal xy tight
xlabel('y')
ylabel('x')
colorbar

