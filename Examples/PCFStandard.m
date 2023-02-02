%% This script computes the modes of a slightly multimode photonic crystal fiber
clear
close all
format long g
[ParentPath, ~, ~] = fileparts(pwd);
addpath([ParentPath '/RefractiveIndexes']);
addpath([ParentPath '/Tools']);

%% Define the simulation window
xmax = 50e-6;
ymax = 50e-6;
N = 512;
x = linspace(-xmax/2, xmax/2, N);
y = x;
[X,Y] = meshgrid(x,y);

%% Define the fiber parameters
lambda = 1064e-9;
fiberParams.Pitch = 3.2e-6;
fiberParams.dop = 1.45/3.2;
fiberParams.Dclad = 35e-6;
fiberParams.lambda = lambda;
RIndexMap = PCFIndex(X, Y, fiberParams);
% Plot RIMap to check if the simulation window is big enough
figure()
imagesc(x*1e6, y*1e6, real(RIndexMap)) % Index is complex since it takes material losses into account
xlabel('x (�m)')
ylabel('y (�m)')
colormap gray
drawnow
%% Run the simulation
nModes = 30; 
n_target = SilicaIndex(lambda); % Since the mode is propagating in fused silica, its effective index can't be higher than the refractive index

[neff, LP] = ModeSolver(RIndexMap, x, y, 'lambda', lambda, 'nModes', ...
    nModes, 'coreRadius', fiberParams.Pitch*2, 'target', n_target, 'plot', true, ...
    'IndexContour', true);