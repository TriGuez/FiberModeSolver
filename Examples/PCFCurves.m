%% This script computes the modes of a slightly multimode photonic crystal fiber
clear
close all
format long g
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Define the simulation window
xmax = 50e-6;
ymax = 50e-6;
N = 512;
x = linspace(-xmax/2, xmax/2, N);
y = x;
[X,Y] = meshgrid(x,y);

%% Define the fiber parameters
pitch = 3.2e-6;
d_over_pitch = 1.4/3.2;
Dclad = 35e-6;
lambda = 1064e-9;
fiberParams = {pitch; d_over_pitch; lambda; Dclad};
RIndexMap = PCFIndex(X, Y, fiberParams);
% Plot RIMap to check if the simulation window is big enough
figure()
imagesc(x*1e6, y*1e6, real(RIndexMap)) % Index is complex since it takes material losses into account
xlabel('x (µm)')
ylabel('y (µm)')
colormap gray

%% Run the simulation
nModes = 1; 
neffs = zeros(1,256);
losses = zeros(1,256);
aeff = zeros(1,256);
lambda = linspace(500e-9,1800e-9,256);
for jk = 1:256
    n_target = SilicaIndex(lambda(jk)); % Since the mode is propagating in fused silica, its effective index can't be higher than the refractive index
    fiberParams = {pitch; d_over_pitch; lambda(jk); Dclad};
    RIndexMap = PCFIndex(X, Y, fiberParams);
    [neff, LP] = ModeSolver(RIndexMap, x, y, 'lambda', lambda(jk), 'nModes', ...
        nModes, 'coreRadius', pitch*2, 'target', n_target, 'plot', false, ...
        'IndexContour', true);
    neffs(jk) = real(neff);
    losses(jk) = (20*imag(neff))./log(10) * 2*pi/lambda(jk);
    aeff(jk) = ModeArea(x, y, LP);
end