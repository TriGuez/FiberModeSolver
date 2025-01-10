function [index] = ZBLANIndex(lbd)

B1 = 1.079256; B2 = 0.152962; B3 = 0.754995;
C1 = 0.008625; C2 = 0.001; C3 = 200; 

lbd = lbd*1e6;

index = sqrt(1+((B1.*lbd.^2)./(lbd.^2-C1))+((B2.*lbd.^2)./(lbd.^2-C2))+((B3.*lbd.^2)./(lbd.^2-C3)));