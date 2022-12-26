function n = StepIndex(X, Y, FiberParams)

%% Description : 
%% This function creates the refractive index profile of step index fiber, including material optical losses (i.e complex refractive index)
%% Inputs : 
%% * X, Y : spatial mesh
%% * FiberParams : Matlab structure with fields : {1} Core diameter (m), {2} Cladding diameter (m), {3} Core numerical apperture , {4} Working wavelength (m)
%%
%% Output : 
%%* n : Refractive index map 

Dcore = FiberParams{1}; Dclad = FiberParams{2}; NA = FiberParams{3};
lambda = FiberParams{4};
alpha = SilicaLosses(lambda)*4.343;
kappa = (alpha.*log(10)/20)*(lambda)/(2*pi);
nClad = SilicaIndex(lambda);
nCore = sqrt(nClad^2 + NA^2);

n = ones(size(X));
n(X.^2+Y.^2<(Dclad/2).^2) = nClad;
n(X.^2+Y.^2<(Dcore/2).^2) = nCore;
n = n + 1i.*kappa;