clc; clear;

nx    = 512;
ny    = 512;
nspat = nx*ny;
nt    = 9;
X     = zeros(nspat, nt);

file_fmt = '%s%05i.bin';
var      = 'density';

for i = 0 : nt-1
    file = sprintf('%s%05i.bin',var, i);
    fid = fopen(file);
    var_data = fread(fid, nspat, 'float64');
    X(:,i+1) = var_data;
    fclose(fid);
end

Phi_py_pre   = load('Phi_PreRsh.txt');
Phi_py       = load('Phi_Compare.txt');
Lambda_py    = load('Lambda_Compare.txt');
Psi_py       = load('Psi_Compare.txt');
R_py         = load('R_Compare.txt');

% ----- Compute POD
R            = X'*X;
[Psi,Lambda] = eig(R);
Phi          = X*Psi/sqrt(Lambda);
A            = X'*Phi;

