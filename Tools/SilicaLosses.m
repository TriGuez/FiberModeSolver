function alpha = SilicaLosses(lbd)
%  This function computes the silica optical losses vs wavelength, according to:
%  Sorensen, S. T. " Deep-blue supercontinuum light sources based on tapered
%  photonic crystal fibers", PhD thesis, DTU Fotonik ,2013.
%  INPUTS :
%        lbd : Wavelength vector [m]
%  OUTPUTS : 
%        losses : Losses vector [m?¹]

auv = 0.001*exp((4.67e-6)./lbd)./1000./4.343;
air = 6e11*exp(-(47.8e-6)./lbd)./1000./4.343;
asc = ((1.3./((lbd*1e6).^4))+1)./1000./4.343;
aoh = (7./(1+((lbd-1.380e-6)./16e-9).^2))./1000./4.343;
alpha = air + auv + asc + aoh;
alpha(lbd<300e-9) = 0;
alpha (lbd > 2500e-9) = 0;