% a example code of 2d two-component Bose-Hubbard model.
%
%
% Notes:
%      |<-       Lx       ->|
%  --- 01  02  03  04  05  06
%   |  07  08  09  10  11  12
%  Ly  13        ...       18     
%   |  19        ...       24
%  --- 25        ...       30
%
%


close all
clc

addpath(genpath('libs'))


%% initial state 
L = 26;
Lx = 13;
Ly = 2;
st_init_up = zeros(Ly,Lx);
st_init_dn = [0 0 1 0 1 0 1 0 1 0 1 0 0;...
    0 0 0 0 0 0 0 0 0 0 0 0 0];

L = 14;
Lx = 14;
Ly = 1;
st_init_up = zeros(Ly,Lx);
st_init_dn = [0 1 1 1 1 1 1 1 1 1 1 1 1 0];

% L = 17;
% Lx = 17;
% Ly = 1;
% st_init_up = zeros(Ly,Lx);
% st_init_dn = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 0 0 0];

if (sum(size(st_init_up) == [Ly,Lx]) ~= 2) ...
        || (sum(size(st_init_dn) == [Ly,Lx]) ~= 2)
    error('Error! Invalid input parameters!')
end

%
N_up = sum(st_init_up(:));
N_dn = sum(st_init_dn(:));
% the maximum occupation number of up/dn state in each lattice site
nMax = 3; 


%% Hubbard parameters
J = 70.0 * 2 * pi;
U = 820.0 * 2 * pi;

%
U_uu = U;      % U_up_up
U_dd = U;      % U_dn_dn
U_ud = U;      % U_up_dn
Jx_up = J * ones(1,Lx-1);
Jx_dn = J * ones(1,Lx-1);
Jy_up = 0;
Jy_dn = 0;

% staggered potential
staG0 = 0.5 * U + 1.0 * J;

mu_up = 0.0*2*pi*(1:Lx) + mod(1:Lx,2)*staG0;
mu_dn = 110.0*2*pi*(1:Lx) + mod(1:Lx,2)*staG0;

mu_up = repmat(mu_up,Ly,1);
mu_dn = repmat(mu_dn,Ly,1);

% boundary conditions
% 'obc' - open boundary condition
% 'pbc' - periodic case
BDCx = 'obc';   % x-direction
BDCy = 'obc';   % y-direction
% BDCx = 'pbc';   % x-direction
% BDCy = 'pbc';   % y-direction


%% evolve time sequence
nt = 101;
T = 100 * 1e-03;          % 100 ms
tl = linspace(0,T,nt);
dt = tl(2) - tl(1);


%% basis generate
fprintf('generating basis.\n')
tBase = tic;
basis = boson_basis_2d_spinor(L,Lx,Ly,N_up,N_dn,nMax);
tD_base = toc(tBase);
fprintf('elapsed time is %.6f seconds.\n',tD_base)
ns = basis.n_bs;
fprintf('total basis number is %d.\n',ns)


%% state index search
idx_init_st = state_index_2d_indicate(st_init_up,st_init_dn,basis);
% basis.state_up(idx_init_st,:)
% basis.state_dn(idx_init_st,:)
phi_init = zeros(ns,1);
phi_init(idx_init_st) = 1;


%% hamiltonian elements generate
fprintf('\ngenerating hamiltonian elements.\n')
tStart = tic; 
ham_elems = hamiltonian_2d_bhm2c_elements(basis,BDCx,BDCy);
tEnd = toc(tStart);
fprintf('total elapsed time is %.6f seconds.\n',tEnd)

fprintf('\ngenerating the hamiltonian.\n')
tStart = tic; 
ham_cur = hamiltonian_2d_bhm2c(basis,ham_elems,Jx_up,Jx_dn,...
    Jy_up,Jy_dn,U_uu,U_dd,U_ud,mu_up,mu_dn);
tEnd = toc(tStart);
fprintf('elapsed time is %.6f seconds.\n',tEnd)
ham = ham_cur.ham;
fprintf('number of nonzero hamiltonian elements: %d.\n',nnz(ham))


%% ground state
k = 3;
[V,D] = eigs(ham,k,'sr');   % find k smallest real eigenvectors and eigenvalues

psi_gs = V(:,1);
energy_gs = D(1,1);

probl_gs = psi_gs .* conj(psi_gs);
figure('Color','w')
plot(probl_gs)
xlim([1,ns])
ax = gca;
ax.FontSize = 14;
xlabel('Basis number','FontSize',16)
ylabel('Probabilility','FontSize',16)


%% evolution
density_up_Mt = [];
density_dn_Mt = [];

psic = phi_init;
probl = abs(conj(psic).*psic);

tic
stat_nC = occupation_statistics_1d_bhm2C_Fcn(psic,basis);
toc
density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);

%
tStart = tic; 
for kk = 2:nt
    fprintf('Current process: %04d / %04d.\n',kk,nt)
    
    tSC = tic;
    psic = expv(-1i*dt,ham,psic,1e-7,30);
    tEC = toc(tSC);
    fprintf('time for evolution: %.6f seconds.\n',tEC)
    
    stat_nC = occupation_statistics_1d_bhm2C_Fcn(psic,basis);
    density_up_Mt = cat(1,density_up_Mt,stat_nC.density_up_ltC);
    density_dn_Mt = cat(1,density_dn_Mt,stat_nC.density_dn_ltC);
    
    % fprintf('\n')
end
tEnd = toc(tStart);
fprintf('\nElapsed time is %.6f seconds.\n',tEnd)


%%
x = 1:L;
y = tl*1000;

figure('Color','w')
imagesc(x,y,density_up_Mt)

figure('Color','w')
imagesc(x,y,density_dn_Mt)



