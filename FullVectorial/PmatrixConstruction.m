function P = PmatrixConstruction(h, Ux, Uy, Vx, Vy, eps_rx, eps_ry, eps_rz,...
    k0)
N = length(Ux(:,1));

Ux = 1/h * Ux;
Uy = 1/h * Uy;
Vx = 1/h * Vx;
Vy = 1/h * Vy;

eps_rzm1 = 1./eps_rz(:);
eps_rzm1(isinf(eps_rzm1)) = 0;

eps_rx = spdiags(eps_rx(:),0,N,N);
eps_ry = spdiags(eps_ry(:),0,N,N);

eps_rzm1 = spdiags(eps_rzm1(:),0,N,N);

Pxx = -(1/(k0*k0))*Ux*eps_rzm1*Vy*Vx*Uy + (k0*k0*speye(N) + ...
    Ux*eps_rzm1*Vx)*(eps_rx + (1/(k0*k0))*Vy*Uy);

Pyy = -(1/(k0*k0))*Uy*eps_rzm1*Vx*Vy*Ux + (k0*k0*speye(N) + ...
    Uy*eps_rzm1*Vy)*(eps_ry+ (1/(k0*k0))*Vx*Ux);

Pxy = Ux*eps_rzm1*Vy*(eps_ry + (1/(k0*k0))*Vx*Ux) - (1/(k0*k0))* ...
    (k0*k0*speye(N) + Ux*eps_rzm1*Vx)*Vy*Ux;

Pyx = Uy*eps_rzm1*Vx*(eps_rx + (1/(k0*k0))*Vy*Uy) - (1/(k0*k0))* ...
    (k0*k0*speye(N) + Uy*eps_rzm1*Vy)*Vx*Uy;

P = [Pxx Pxy; Pyx Pyy];