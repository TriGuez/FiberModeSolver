function [eps_rx, eps_ry, eps_rz] = RIndexMoy(eps_r)

N = length(eps_r(1,:));
eps_rx = (eps_r + circshift(eps_r,-1,2))./2;
eps_ry = (eps_r + circshift(eps_r,-1,1))./2;
eps_rz = (eps_r + circshift(eps_r,[-1 -1]) + circshift(eps_r,-1,2) + ...
     circshift(eps_r,-1,1))./4;
% eps_rx = spdiags(((eps_rx(:)).^2),0,N^2,N^2);
% eps_ry = spdiags(((eps_ry(:)).^2),0,N^2,N^2);
% eps_rz = spdiags(((eps_rz(:)).^2),0,N^2,N^2);