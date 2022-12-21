function n = ARFIndex(X,Y, nParameters)

Dclad = nParameters{1};
Ncap = nParameters{2};
DextCap = nParameters{3};
Ecap = nParameters{4};
lbd = nParameters{5};

Rclad = Dclad/2;
Rcap = DextCap/2;
RintCap = Rcap - Ecap;
alpha = silicaLosses(lbd)*4.343;
% kappa = (alpha.*log(10)/20)*(nParameters{5})/(2*pi);
kappa = 0;
n = (silica_index(lbd)+1i*kappa).*ones(size(X));
n(X.^2+Y.^2<(Dclad/2).^2) = 1;
alpha = 2*pi/Ncap;
for jk = 1:Ncap
    n((X-(Rclad-Rcap).*cos(jk*alpha)).^2 + (Y-(Rclad-Rcap).*sin(jk*alpha)).^2 < Rcap.^2) = silica_index(lbd)+1i*kappa;
    n((X-(Rclad-Rcap).*cos(jk*alpha)).^2 + (Y-(Rclad-Rcap).*sin(jk*alpha)).^2 < RintCap.^2) = 1;
end
n(X.^2+Y.^2 > (Rclad+20e-6).^2) = 1;