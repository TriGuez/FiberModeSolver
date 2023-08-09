function n = PCFIndex(X,Y, FiberParams)
%% This function creates the refractive index profile of an air silica photonic crystal fiber (One hole missing triangular lattice, 5 rows of holes), including material optical losses (i.e complex refractive index)
%% Inputs : 
%% * X, Y : spatial mesh
%% * FiberParams : Matlab structure with fields : {1} lattice pitch (m), {2} $d\over \Lambda$ parameter, {3} Working wavelength (m), {4} Cladding diameter
%%
%% Output : 
%% * n : Refractive index map 

format long g
pitch = FiberParams.Pitch;
d_over_pitch = FiberParams.dop;
Dclad = FiberParams.Dclad;
n_anneaux = FiberParams.n_anneaux;
d = d_over_pitch.*pitch;
lambda = FiberParams.lambda;
n_background = SilicaIndex(lambda);
n = zeros(size(X));
n(X.^2+Y.^2<(Dclad/2).^2) = n_background;
alpha = SilicaLosses(lambda)*4.343;
kappa = (alpha.*log(10)/20)*(lambda)/(2*pi);
n = n + 1i*kappa;


x_holes = zeros(1,6*n_anneaux*(n_anneaux+1)/2);
y_holes = zeros(1,6*n_anneaux*(n_anneaux+1)/2);
nbre_cyl = zeros(1,n_anneaux);


for indice_anneau = 1 : n_anneaux

    nbre_cyl(indice_anneau) = 0;
    T = zeros(2,4*indice_anneau+1);
    for m = 1 : indice_anneau
        T(1,m+1)     = -0.5;
        T(1,m+1+indice_anneau)   = -1;
        T(1,m+1+2*indice_anneau) = -0.5;
    end
    for m = 1 : 3*indice_anneau
        T(1,m+1+3*indice_anneau) = -T(1,m+1);
    end
    for m = 1 : indice_anneau
        T(2,m+1)     = sqrt(3)/2;
        T(2,m+1+indice_anneau)   = 0;
        T(2,m+1+2*indice_anneau) = -sqrt(3)/2;
    end
      
    for j = 1 : 3*m
        T(2,m+1+3*indice_anneau) = -T(2,m+1);
    end
       
    coord_x = 0;
    coord_y = 0;
    for k = 1 : length(T(1,:))/3
        coord_x = T(1,k) + coord_x;
        coord_y = coord_y + T(2,k);
        if ~( (m + coord_x) == 0)
            theta = atan(coord_y/(m + coord_x));
        else
            theta = pi/2;
        end
        if (coord_x == 0)
            theta = 0;
            rho = m;
        else
            rho = coord_y/sin(theta);
        end

        if  ( rho < 0 )
            rho = abs(rho);
            theta = pi + theta;
        end
            
        for il = 1 : 3
            nbre_cyl(indice_anneau) = nbre_cyl(indice_anneau) + 1;
            x_holes(sum(nbre_cyl)) = rho*pitch*cos(theta);
            y_holes(sum(nbre_cyl)) = rho*pitch*sin(theta);
            theta = theta + (2*pi/3);
        end
       
    end
end

for jk = 1:length(x_holes)
    n((X-x_holes(jk)).^2 + (Y-y_holes(jk)).^2 < (d/2)^2) = 1;
end











% % first ring
% alpha = pi/3;
% for jk = 1:6
%     n((X-(pitch.*cos(jk*alpha))).^2 + (Y-(pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
% end
% 
% % 2nd ring
% alpha = pi/6;
% for jk = 1:12
%     if mod(jk,2)==0
%         n((X-(2*pitch.*cos(jk*alpha))).^2 + (Y-(2*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
%     else
%         n((X-(sqrt(3)*pitch.*cos(jk*alpha))).^2 + (Y-(sqrt(3)*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
%     end
% end
% 
% % 3rd ring
% alpha = pi/3;
% 
% for jk = 1:6
%     n((X-(3*pitch.*cos(jk*alpha))).^2 + (Y-(3*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
% end
% 
% n((X-(2.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(2*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% 
% n((X+(2.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X+(2*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% 
% n((X-(2.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(2*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% 
% n((X+(2.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X+(2*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% 
% n((X-(pitch.*cos(5*pi/3))).^2 + (Y-(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(4*pi/3))).^2 + (Y-(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% 
% n((X-(pitch.*cos(5*pi/3))).^2 + (Y+(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(4*pi/3))).^2 + (Y+(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% 
% % 4th ring
% 
% for jk = 1:6
%     n((X-(4*pitch.*cos(jk*alpha))).^2 + (Y-(4*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
% end
% 
% n((X-(3.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X-(2.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% 
% n((X+(3.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X+(3*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X+(2.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% 
% n((X+(3.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X+(3*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X+(2.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% 
% n((X-(3.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X-(2.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% 
% 
% n((X-(pitch.*cos(6*pi/3))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(3*pi/3))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(pi/2))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% 
% 
% n((X-(pitch.*cos(6*pi/3))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(3*pi/3))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(pi/2))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% 
% 
% % 5th ring
% for jk = 1:6
%     n((X-(5*pitch.*cos(jk*alpha))).^2 + (Y-(5*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
% end
% 
% n((X-(4.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(4*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X-(3.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch)).^2 + (Y-(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;
% 
% n((X+(4.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X+(4*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X+(3.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% n((X+(3*pitch)).^2 + (Y-(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;
% 
% n((X-(4.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(4*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X-(3.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch)).^2 + (Y+(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;
% 
% n((X+(4.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
% n((X+(4*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
% n((X+(3.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
% n((X+(3*pitch)).^2 + (Y+(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;
% 
% n((X-(pitch.*cos(5*pi/3))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(2*pi/3))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch.*cos(2*pi/6))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch.*cos(4*pi/6))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% 
% n((X-(pitch.*cos(5*pi/3))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(pitch.*cos(2*pi/3))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch.*cos(2*pi/6))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
% n((X-(3*pitch.*cos(4*pi/6))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;

