%% This script computes the modes of a AntiResonnant Fiber
clear
close all
format long g
[ParentPath, ~, ~] = fileparts(pwd);
addpath([ParentPath '/RefractiveIndexes']);
addpath([ParentPath '/Tools']);


%% Define the simulation window
xmax = 150e-6;
ymax = 150e-6;
N = 512;
x = linspace(-xmax/2, xmax/2, N);
y = x;
[X,Y] = meshgrid(x,y);

%% Define the fiber parameters
lambda = 1550e-9;
fiberParams.lambda = lambda;
fiberParams.Dcore = 105e-6;
fiberParams.Dclad = 125e-6;
fiberParams.NA = 0.22;

RIndexMap = StepIndex(X, Y, fiberParams);
% Plot RIMap to check if the simulation window is big enough
figure()
imagesc(x*1e6, y*1e6, real(RIndexMap)) % Index is complex since it takes material losses into account
xlabel('x (µm)')
ylabel('y (µm)')
colormap gray

%% Run the simulation
nModes = 5; 
n_target = sqrt(SilicaIndex(lambda)^2+fiberParams.NA^2); % Since the mode is propagating in the core, its effective index can't be higher than the refractive index

[neff, LP] = ModeSolver(RIndexMap, x, y, 'lambda', lambda, 'nModes', ...
    nModes, 'coreRadius', fiberParams.Dcore/2, 'target', n_target, 'plot', true, ...
    'IndexContour', true);