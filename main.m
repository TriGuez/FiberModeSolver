clear
close all
format long g
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

tic
xmax = 150e-6;
ymax = 150e-6;
N = 512;
h = (xmax)/N;
x = linspace(-xmax/2, xmax/2, N);
y = x;
[X,Y] = meshgrid(x,y);
% lbd = linspace(500e-9,1800e-9,256);
lbd = 1064e-9;
neff_sim = zeros(1,length(lbd));
Aeff = zeros(1,length(lbd));
nvis = zeros(1,length(lbd));
MFD = zeros(1,length(lbd));
LP = zeros(N,N,length(lbd));
for jk = 1:length(lbd)

    lambda = lbd(jk);
    
%     nParameters = {3.2e-6;1.4/3.2; lambda; 35e-6};
%     nParameters = {97e-6; 6; 26e-6; 1150e-9; lambda}; %ARF
    nParameters = {105e-6; 125e-6; 0.22; lambda};
    n = StepIndex(X,Y,nParameters);
%     n = PCFIndex(X,Y,nParameters);
%     n = ARFIndex(X,Y,nParameters);
    
    figure()
    imagesc(real(n));
    colormap gray
    
    [neffs, LP] = ModeSolver(n,x,y,'lambda', lambda,...
        'nModes', 50, 'coreRadius', 55e-6, 'target',sqrt(SilicaIndex(lambda)^2+0.22^2) ...
        , 'plot',true, 'IndexContour', true);
%     neff_sim(jk) = neffs;
%     Aeff(jk) = ModeArea(x,y,LP(:,:,jk));
%     MFD(jk) = ModeFieldDiameter(x,y,LP(:,:,jk),'4sigma');
end
toc