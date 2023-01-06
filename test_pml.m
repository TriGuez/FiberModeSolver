%% This script computes the modes of a slightly multimode photonic crystal fiber
clear
close all
format long g
[ParentPath, ~, ~] = fileparts(pwd);
addpath([ParentPath '/RefractiveIndexes']);
addpath([ParentPath '/Tools']);

%% Define the simulation window
xmax = 35e-6;
ymax = 35e-6;
N = 256;
% x = linspace(-xmax/2, xmax/2, N);
pad = 2.0;
Nt = round((pad*xmax)/(xmax/N));
x = linspace(-(xmax/2)*pad,(xmax/2)*pad,Nt);
y = x;
[X,Y] = meshgrid(x,y);

%% Define the fiber parameters
lambda = 500e-9;
fiberParams.Pitch = 3.2e-6;
fiberParams.dop = 1.45/3.2;
fiberParams.Dclad = 35e-6;
% fiberParams.Dcore = 40e-6;
% fiberParams.NA = 0.045;
% fiberParams.Dclad = 97e-6;
% fiberParams.NCap = 6;
% fiberParams.DextCap = 26e-6;
% fiberParams.ECap = 1150e-9;
fiberParams.lambda = lambda;
% RIndexMap = ARFIndex(X, Y, fiberParams);
RIndexMap = PCFIndex(X, Y, fiberParams);
% RIndexMap = StepIndex(X, Y, fiberParams);
BendingParams.R = 150e-3;
BendingParams.Angle = 0;
% RIndexMap = BendedIndex(X, Y, RIndexMap, BendingParams);
% Plot RIMap to check if the simulation window is big enough
figure()
imagesc(x*1e6, y*1e6, real(RIndexMap)) % Index is complex since it takes material losses into account
xlabel('x (�m)')
ylabel('y (�m)')
colormap gray

%% Define the PML
gamma = 1i-1;
sx = 1-gamma.*((3*lambda)./(4*pi*(xmax/2))).*(X./(xmax/2)).^2.*log(1e-12);
sx(abs(X)<xmax/2) = 1;
sy = 1-gamma.*((3*lambda)./(4*pi*(xmax/2))).*(Y./(xmax/2)).^2.*log(1e-12);
sy(abs(Y)<xmax/2) = 1;
% sx= ones(size(X));
% sy=sx;
RIndexMap = sx.*RIndexMap;
RIndexMap = sy.*RIndexMap;

%% Run the simulation
nModes = 10; 
n_target = ieim_pcf(fiberParams.Pitch, fiberParams.dop,fiberParams.lambda,SilicaIndex(lambda));
% n_target = 1;
% n_target = sqrt(SilicaIndex(lambda).^2 + fiberParams.NA.^2); % Since the mode is propagating in fused silica, its effective index can't be higher than the refractive index

[neff, LP] = ModeSolver(RIndexMap, x, y, 'lambda', lambda, 'nModes', ...
    nModes, 'coreRadius',(35/2)*1e-6, 'target', n_target, 'plot', true, ...
    'IndexContour', true, 'pml', cat(3,sx,sy));
losses = 20.*imag(neff)./log(10) .* (2*pi./lambda) * 1000