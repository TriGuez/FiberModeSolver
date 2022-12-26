function n = ARFIndex(X,Y, FiberParams)
%% This function creates the refractive index profile of a fused silica made anti-resonnant fiber, including material optical losses (i.e complex refractive index)
%% Inputs : 
%% - X, Y : spatial grid
%% - FiberParams : Matlab struct with fields : {1} Cladding diameter, {2} number of capillaries, {3} Outter capillaries diameter, {4} capillaries thickness, {5} Working wavelength
%% Outputs : 
%% - n : Refractive index profile of the fiber

Dclad = FiberParams{1};
Ncap = FiberParams{2};
DextCap = FiberParams{3};
Ecap = FiberParams{4};
lbd = FiberParams{5};

Rclad = Dclad/2;
Rcap = DextCap/2;
RintCap = Rcap - Ecap;
alpha = SilicaLosses(lbd)*4.343;
kappa = (alpha.*log(10)/20)*(FiberParams{5})/(2*pi);
n = (SilicaIndex(lbd)+1i*kappa).*ones(size(X));
n(X.^2+Y.^2<(Dclad/2).^2) = 1;
alpha = 2*pi/Ncap;

for jk = 1:Ncap
    n((X-(Rclad-Rcap).*cos(jk*alpha)).^2 + (Y-(Rclad-Rcap).*sin(jk*alpha)).^2 < Rcap.^2) = SilicaIndex(lbd)+1i*kappa;
    n((X-(Rclad-Rcap).*cos(jk*alpha)).^2 + (Y-(Rclad-Rcap).*sin(jk*alpha)).^2 < RintCap.^2) = 1;
end
n(X.^2+Y.^2 > (Rclad+20e-6).^2) = 1;