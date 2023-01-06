function n = BendedIndex(X ,Y , RIndexMap, BendingParams)


%% This function computes the impact of bending on a given refractive index profile according to [3].
%%
%% Inputs : 
%% * X, Y : spatial mesh
%% * RIndexMap : Input refractive index profile
%% * BendingParams : Matlab structure with fields : R : Bending radius (m), Angle : Angle between bending direction & the horizontal direction (rad).
%% 
%% Outputs : 
%% * n : Bended refractive index profile

rho_e = 0.22; % Effect of photoelastic tensor & Poisson's ratio in fused silica
R = BendingParams.R;
Angle = BendingParams.Angle;
x = X(1,:);
y = Y(:,1);
bendDir = x.*cos(Angle) - y.*sin(Angle);
kappa = imag(RIndexMap);
RIndexMap = real(RIndexMap);

BendMap = RIndexMap.*(1-RIndexMap.^2.*(bendDir./(2*R)).*rho_e).*exp(bendDir/R);

n = ones(size(X)) + 0*1i;
n((RIndexMap) > 1) = BendMap((RIndexMap) > 1);
n = n + 1i.*kappa;
