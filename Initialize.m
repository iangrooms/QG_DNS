% Set simulation parameters
N = 256; % Number of points in each direction
countDiag = 5; % Compute diagnostics every countDiag steps
isAdaptive = false; % Is the time step adaptive
dt = .01; % time step size; if adaptive then it will change
Nt = 1E3; % Number of time steps
qlim = 1E3; % if any q > qlim, simulation stops

% Set physical parameters
kd = 1; % Nondimensional deformation wavenumber
LX = 2*pi*16; % Nondimensional domain width
kb = 0; % Nondimensional beta wavenumber
r = 0; % Nondimensional Ekman friction coefficient
c_d = 0.1; % Nondimensional quadratic drag coefficient
nu = 5E-7; % Coefficient of hyperviscous vorticity diffusion

% Put useful stuff into a struct
params = struct('kd',kd,'kb',kb,'r',r,'cd',c_d,'nu',nu,'N',N,'LX',LX);

% Set up hyperviscous PV dissipation
k = (2*pi/LX)*[0:N/2 -N/2+1:-1]'; % wavenumbers
L = zeros([N N 2]);
for jj=1:N
    for ii=1:N
        kr = sqrt(k(ii)^2+k(jj)^2);
        L(ii,jj,:) = -nu*kr^8;
    end
end
clear k kr ii jj lf

% adaptive stepping stuff:
tol= 1E-2;
r0 = .8*tol;

% Initialize
t = 0;
% qp is the "physical" potential vorticity, i.e. values on the grid
qp = randn([N N 2]); 
% q is the Fourier coefficients of potential vorticity
q = fft2(qp);

% Initialize diagnostics
T = zeros(1,floor(Nt/countDiag));
