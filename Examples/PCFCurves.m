%% This script gives an example of a multi-wavelength simulation
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
fiberParams.dop = 1.4/3.2;
fiberParams.Dclad = 35e-6;
fiberParams.lambda = 1064e-9;
RIndexMap = PCFIndex(X, Y, fiberParams);
% Plot RIMap to check if the simulation window is big enough
figure()
imagesc(x*1e6, y*1e6, real(RIndexMap)) % Index is complex since it takes material losses into account
xlabel('x (�m)')
ylabel('y (�m)')
colormap gray

%% Run the simulation
nModes = 1; 
neffs = zeros(1,256);
losses = zeros(1,256);
aeff = zeros(1,256);
lambda = linspace(500e-9,1800e-9,256);
for jk = 1:256
    n_target = SilicaIndex(lambda(jk)); % Since the mode is propagating in fused silica, its effective index can't be higher than the refractive index
    fiberParams.lambda = lambda(jk);
    RIndexMap = PCFIndex(X, Y, fiberParams);
    [neff, LP] = ModeSolver(RIndexMap, x, y, 'lambda', lambda(jk), 'nModes', ...
        nModes, 'coreRadius', fiberParams.Pitch*2, 'target', n_target, 'plot', false, ...
        'IndexContour', true);
    neffs(jk) = real(neff);
    losses(jk) = (20*imag(neff))./log(10) * 2*pi/lambda(jk);
    aeff(jk) = ModeArea(x, y, LP);
end

%% Plotting results 
figure()
plot(lambda.*1e9, neff)
xlabel('Wavelength (nm)')
ylabel('n_{eff}')
grid on
grid minor
figure()
yyaxis left
plot(lambda*1e9, aeff.*1e12)
xlabel('Wavelength')
ylabel('A_{eff} (µm^{2})')
yyaxis right
plot(lambda*1e9, losses.*4.343)
ylabel('Attenuation (dB.m^{-1})')