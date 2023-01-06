function n = ParabolicIndex(X, Y, FiberParams)
    
%% This function creates the refractive index profile of a parabolic graded index fiber, including material losses (i.e complex refractive index)
%% 
%% Inputs : 
%% * X, Y : spatial mesh
%% * FiberParams : Matlab structure with fields : Dcore : Core diameter (m), Dclad : Cladding diameter (m), dn : Refractive index difference of the parabola, lambda : WOrking wavelength (m)
%% 
%% Outputs : 
%% * n : refractive index profile
   
    CoreRadius = FiberParams.Dcore/2;
    dn = FiberParams.dn;
    lambda = FiberParams.lambda;
    Dclad = FiberParams.Dclad;
    alpha = SilicaLosses(lbd)*4.343;
    kappa = (alpha.*log(10)/20)*(lbd)/(2*pi);

    n = ones(size(X));
    n(X.^2+Y.^2<(Dclad/2).^2) = SilicaIndex(lambda) + 1i.*kappa;
    R = sqrt((X).^2 + (Y).^2);
    n(R < CoreRadius) = SilicaIndex(lambda) + (dn)*(1 - (R(R < CoreRadius)/CoreRadius).^2) + 1i.*kappa; % Equation for parabolic graded index core
end