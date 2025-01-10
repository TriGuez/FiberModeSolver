%% This script computes the modes of a step_index ZBLAN fiber
clear
close all
format long g
[ParentPath, ~, ~] = fileparts(pwd);
addpath([ParentPath '/RefractiveIndexes']);
addpath([ParentPath '/Tools']);


%% Define the simulation window
xmax = 150e-6;
ymax = 150e-6;
N = 256;
x = linspace(-xmax/2, xmax/2, N);
y = x;
[X,Y] = meshgrid(x,y);

%% Define the fiber parameters
%lbd = linspace(2000,4000, 256)*1e-9;
lbd = 3000e-9;
%neff = zeros(1,256);
for jk = 1:length(lbd)
fiberParams.lambda = lbd(jk);
fiberParams.Dcore = 15e-6;
fiberParams.Dclad = 125e-6;
fiberParams.NA = 0.15;

RIndexMap = StepIndexZBLAN(X, Y, fiberParams);
% Plot RIMap to check if the simulation window is big enough
% figure()
% imagesc(x*1e6, y*1e6, real(RIndexMap)) 
% xlabel('x (µm)')
% ylabel('y (µm)')
% colormap gray

%% Run the simulation
nModes = 1; 
n_target = sqrt(ZBLANIndex(lbd(jk))^2+fiberParams.NA^2); % Since the mode is propagating in the core, its effective index can't be higher than the refractive index

[neff(jk), LP] = ModeSolver(RIndexMap, x, y, 'lambda', lbd(jk), 'nModes', ...
    nModes, 'coreRadius', fiberParams.Dcore/2, 'target', n_target, 'plot', true, ...
    'IndexContour', true);
end