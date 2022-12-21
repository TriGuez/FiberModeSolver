%% This script computes the modes of a AntiResonnant Fiber
clear
close all
format long g
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

%% Define the simulation window
xmax = 150e-6;
ymax = 150e-6;
N = 512;
x = linspace(-xmax/2, xmax/2, N);
y = x;
[X,Y] = meshgrid(x,y);

%% Define the fiber parameters
lambda = 1550e-9;
Dclad = 95e-6;
NCapillaries = 5;
DCapillaries = 25e-6;
ECapillaries = 1.5e-6;
fiberParams = {Dclad; NCapillaries; DCapillaries; ECapillaries; lambda}; %ARF
RIndexMap = ARFIndex(X, Y, fiberParams);
% Plot RIMap to check if the simulation window is big enough
figure()
imagesc(x*1e6, y*1e6, real(RIndexMap)) % Index is complex since it takes material losses into account
xlabel('x (µm)')
ylabel('y (µm)')
colormap gray

%% Run the simulation
nModes = 50; 
n_target = 1; % Since the mode is propagating in air, its effective index can't be higher than the refractive index

[neff, LP] = ModeSolver(RIndexMap, x, y, 'lambda', lambda, 'nModes', ...
    nModes, 'coreRadius', 50e-6, 'target', n_target, 'plot', true, ...
    'IndexContour', true);