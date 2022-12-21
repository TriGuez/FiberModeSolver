function [index] = SilicaIndex(lambda) 

lambda = lambda*1e6;

sl1= 68.4043*1e-3;  sa1=0.6961663; 
sl2= 116.2414*1e-3; sa2=0.4079426; 
sl3= 9896.161*1e-3; sa3=0.8974794; 

X1 = sl1./lambda; X2 = sl2./lambda; X3 = sl3./lambda; 
X1c = X1.*X1; X2c = X2.*X2; X3c = X3.*X3; 
T1 = 1-X1c; T2 = 1-X2c; T3 = 1-X3c; 
 
silica_epsilon = 1+sa1./T1+sa2./T2+sa3./T3;
index = sqrt(silica_epsilon);
