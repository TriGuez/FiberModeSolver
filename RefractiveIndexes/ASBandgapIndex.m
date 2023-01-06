function n = ASBandgapIndex(X,Y, FiberParams)

%% This function creates the refractive index profile of an all-solid photonic bandgap fiber (hexagonal lattice of high-index parabolic inclusions), including material optical losses (i.e complex refractive index)
%% 
%% Inputs : 
%% * X, Y : spatial mesh
%% * FiberParams : Matlab structure with fields : Pitch : lattice pitch (m), dop : $d\over \Lambda$ parameter, dn : refractive index difference of the high index inclusions, Dclad : Cladding diameter, lambda : working wavelength (m)
%% 
%% Outputs : 
%% * n : Refractive index map

format long g
pitch = FiberParams.Pitch;
d_over_pitch = FiberParams.dop;
Dclad = FiberParams.Dclad;
d = d_over_pitch.*pitch;
n_background = SilicaIndex(FiberParams.lambda);
n = ones(size(X));
n(X.^2+Y.^2<(Dclad/2).^2) = n_background;
lambda = FiberParams.lambda;
alpha = SilicaLosses(lambda)*4.343;
kappa = (alpha.*log(10)/20)*(lambda)/(2*pi);

r = d/2;
dn = FiberParams.dn;
% first ring
alpha = pi/3;
for jk = 1:6
    R = sqrt((X - (pitch.*cos(jk*alpha))).^2 + (Y - (pitch.*sin(jk*alpha))).^2);
    n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
end
% 
% 2nd ring
alpha = pi/6;
for jk = 1:12
    if mod(jk,2)==0
        R = sqrt((X-(2*pitch.*cos(jk*alpha))).^2 + (Y-(2*pitch.*sin(jk*alpha))).^2);
        n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
    else
        R = sqrt((X-(sqrt(3)*pitch.*cos(jk*alpha))).^2 + (Y-(sqrt(3)*pitch.*sin(jk*alpha))).^2);
        n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
    end
end

% 3rd ring
alpha = pi/3;

for jk = 1:6
    R = sqrt((X-(3*pitch.*cos(jk*alpha))).^2 + (Y-(3*pitch.*sin(jk*alpha))).^2);
    n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
end
R = sqrt((X-(2.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(2*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X+(2.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(2*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(2.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(2*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X+(2.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(2*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(pitch.*cos(5*pi/3))).^2 + (Y-(3*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(4*pi/3))).^2 + (Y-(3*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(pitch.*cos(5*pi/3))).^2 + (Y+(3*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(4*pi/3))).^2 + (Y+(3*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

% 4th ring

for jk = 1:6
    R = sqrt((X-(4*pitch.*cos(jk*alpha))).^2 + (Y-(4*pitch.*sin(jk*alpha))).^2);
    n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
end

R = sqrt((X-(3.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(2.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X+(3.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(3*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(2.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X+(3.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(3*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(2.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 );
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(3.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(2.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(pitch.*cos(6*pi/3))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(3*pi/3))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(pi/2))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(pitch.*cos(6*pi/3))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(3*pi/3))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(pi/2))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

% 5th ring
for jk = 1:6
    R = sqrt((X-(5*pitch.*cos(jk*alpha))).^2 + (Y-(5*pitch.*sin(jk*alpha))).^2);
    n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
end

R = sqrt((X-(4.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(4*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch)).^2 + (Y-(4*pitch*sin(-5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X+(4.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(4*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(3.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(3*pitch)).^2 + (Y-(4*pitch*sin(-5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(4.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(4*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(4*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch)).^2 + (Y+(4*pitch*sin(-5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X+(4.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(4*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(3.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(3*pitch)).^2 + (Y+(4*pitch*sin(-5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(pitch.*cos(5*pi/3))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(2*pi/3))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch.*cos(2*pi/6))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch.*cos(4*pi/6))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);

R = sqrt((X-(pitch.*cos(5*pi/3))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(pitch.*cos(2*pi/3))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch.*cos(2*pi/6))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X-(3*pitch.*cos(4*pi/6))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
n = n + 1i*kappa;