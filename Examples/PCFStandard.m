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
d_over_pitch = 1.8/3.2;
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
nModes = 30; 
n_target = SilicaIndex(lambda); % Since the mode is propagating in fused silica, its effective index can't be higher than the refractive index

[neff, LP] = ModeSolver(RIndexMap, x, y, 'lambda', lambda, 'nModes', ...
    nModes, 'coreRadius', pitch*2, 'target', n_target, 'plot', true, ...
    'IndexContour', true);