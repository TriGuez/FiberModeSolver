clear
close all
format long g
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

tic
xmax = 50e-6;
ymax = 50e-6;
N = 512;
h = (xmax)/N;
x = linspace(-xmax/2, xmax/2, N);
y = x;
[X,Y] = meshgrid(x,y);
% lbd = linspace(500e-9,1800e-9,256);
lbd = 1064e-9
neff_sim = zeros(1,length(lbd));
Aeff = zeros(1,length(lbd));
nvis = zeros(1,length(lbd));
MFD = zeros(1,length(lbd));
LP = zeros(N,N,length(lbd));
for jk = 1:length(lbd)

    lambda = lbd(jk);
    
    nParameters = {3.2e-6;1.4/3.2; lambda; 35e-6};
    % nParameters = {97e-6; 6; 26e-6; 1150e-9; lambda}; %ARF
    n = PCFIndex(X,Y,nParameters);
    % n = ARFIndex(X,Y,nParameters);
    
%     figure()
%     imagesc(real(n));
%     colormap gray
    
    [neffs, LP(:,:,jk)] = ModeSolver(n,x,y,'lambda', lambda,...
        'nModes', 1, 'coreRadius', 2.5e-6, 'target',SilicaIndex(lambda)...
        , 'plot',true, 'IndexContour', false);
    neff_sim(jk) = neffs;
    Aeff(jk) = modeArea(x,y,LP(:,:,jk));
    MFD(jk) = ModeFieldDiameter(x,y,LP(:,:,jk),'4sigma');
    nvis(jk) = ieim_pcf(nParameters{1},nParameters{2}, lambda, 'silica');
end
plot(lbd,neff_sim, '-blue', lbd, nvis, '-red')
toc