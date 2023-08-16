function [neffs, Ex, Ey, Ez, Hx, Hy, Hz, S] = PostProcess(V, D, k0, x, y, eps_rx, eps_ry, eps_rz)

h = max(x(2)-x(1), y(2)-y(1));
N = sqrt(length(V(:,1))/2);
[Ux, Uy, Vx, Vy] = UVMatrixConstruction(N);
eps_rx = spdiags(eps_rx(:),0,N^2,N^2);
eps_ry = spdiags(eps_ry(:),0,N^2,N^2);

eps_rzm1 = 1./eps_rz(:);
eps_rzm1(isinf(eps_rzm1)) = 0;
eps_rzm1 = spdiags(eps_rzm1(:),0,N^2,N^2);

Ux = 1/h * Ux;
Uy = 1/h * Uy;
Vx = 1/h * Vx;
Vy = 1/h * Vy;

neffs = sqrt(diag(D)/k0^2);
idxs = find(~isnan(neffs));
neffs = neffs(idxs);
betas = neffs*k0;

Ex = reshape(V(1:N^2,idxs), N, N, length(idxs));
Ey = reshape(V(N^2+1:end,idxs), N, N, length(idxs));

for jk = 1:length(idxs)
    Extemp = Ex(:,:,jk);
    Extemp = spdiags(Extemp(:),0,N^2,N^2);
    Eytemp = Ey(:,:,jk);
    Eytemp = spdiags(Eytemp(:),0,N^2,N^2);
    Hztemp = (-Uy*Extemp + Ux*spdiags(Eytemp(:),0,N^2,N^2))./(1i*k0);
    Hytemp = (1/(1i*betas(jk))) .* (Vy*Hztemp + 1i*k0.*eps_rx*Extemp);
    Hxtemp = (1/(1i*betas(jk))) .* (Vx*Hztemp - 1i*k0.*eps_ry*Eytemp);
    Eztemp = (-1/(1i*k0)).*eps_rzm1*(-Vy*Hxtemp + Vx*Hytemp);
    Sztemp = Extemp*Hytemp-Eytemp*Hxtemp;
    Sytemp = Eztemp*Hxtemp - Extemp*Hztemp;
    Sxtemp = Eytemp*Hztemp - Eztemp*Hytemp;
    
    Hxtemp = diag(Hxtemp);
    Hx(:,:,jk) = reshape(full(Hxtemp),N,N);
    Hytemp = diag(Hytemp);
    Hy(:,:,jk) = reshape(full(Hytemp),N,N);    
    Hztemp = diag(Hztemp);
    Hz(:,:,jk) = reshape(full(Hztemp),N,N);
    Eztemp = diag(Eztemp);
    Ez(:,:,jk) = reshape(full(Eztemp),N,N);
    Sztemp = diag(Sztemp);
    Sz(:,:,jk) = reshape(full(Sztemp),N,N);
    Sytemp = diag(Sytemp);
    Sy(:,:,jk) = reshape(full(Sytemp),N,N);
    Sxtemp = diag(Sxtemp);
    Sx(:,:,jk) = reshape(full(Sxtemp),N,N);
    S(:,:,jk) = sqrt(Sx(:,:,jk).*conj(Sx(:,:,jk)) + Sy(:,:,jk).*conj(Sy(:,:,jk)) + Sz(:,:,jk).*conj(Sz(:,:,jk)));
end