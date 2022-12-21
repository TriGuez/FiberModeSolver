function n = PCFIndex(X,Y, nParameters)
format long g
pitch = nParameters{1};
d_over_pitch = nParameters{2};
dclad = nParameters{4};
d = d_over_pitch.*pitch;
n_background = SilicaIndex(nParameters{3});
n = ones(size(X));
n(X.^2+Y.^2<(dclad/2).^2) = n_background;
alpha = SilicaLosses(nParameters{3})*4.343;
kappa = (alpha.*log(10)/20)*(nParameters{3})/(2*pi);
n = n + 1i*kappa;

% first ring
alpha = pi/3;
for jk = 1:6
    n((X-(pitch.*cos(jk*alpha))).^2 + (Y-(pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
end

% 2nd ring
alpha = pi/6;
for jk = 1:12
    if mod(jk,2)==0
        n((X-(2*pitch.*cos(jk*alpha))).^2 + (Y-(2*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
    else
        n((X-(sqrt(3)*pitch.*cos(jk*alpha))).^2 + (Y-(sqrt(3)*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
    end
end

% 3rd ring
alpha = pi/3;

for jk = 1:6
    n((X-(3*pitch.*cos(jk*alpha))).^2 + (Y-(3*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
end

n((X-(2.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X-(2*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;

n((X+(2.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X+(2*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;

n((X-(2.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X-(2*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;

n((X+(2.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X+(2*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;

n((X-(pitch.*cos(5*pi/3))).^2 + (Y-(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(4*pi/3))).^2 + (Y-(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;

n((X-(pitch.*cos(5*pi/3))).^2 + (Y+(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(4*pi/3))).^2 + (Y+(3*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;

% 4th ring

for jk = 1:6
    n((X-(4*pitch.*cos(jk*alpha))).^2 + (Y-(4*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
end

n((X-(3.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X-(3*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X-(2.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;

n((X+(3.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X+(3*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X+(2.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;

n((X+(3.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X+(3*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X+(2.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;

n((X-(3.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X-(3*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X-(2.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;


n((X-(pitch.*cos(6*pi/3))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(3*pi/3))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(pi/2))).^2 + (Y+(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;


n((X-(pitch.*cos(6*pi/3))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(3*pi/3))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(pi/2))).^2 + (Y-(4*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;


% 5th ring
for jk = 1:6
    n((X-(5*pitch.*cos(jk*alpha))).^2 + (Y-(5*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
end

n((X-(4.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X-(4*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X-(3.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
n((X-(3*pitch)).^2 + (Y-(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;

n((X+(4.5*pitch)).^2 + (Y-(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X+(4*pitch)).^2 + (Y-(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X+(3.5*pitch)).^2 + (Y-(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
n((X+(3*pitch)).^2 + (Y-(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;

n((X-(4.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X-(4*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X-(3.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
n((X-(3*pitch)).^2 + (Y+(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;

n((X+(4.5*pitch)).^2 + (Y+(pitch*sin(2*pi/3))).^2 < (d/2).^2) = 1;
n((X+(4*pitch)).^2 + (Y+(sqrt(3)*pitch*sin(3*pi/6))).^2 < (d/2).^2) = 1;
n((X+(3.5*pitch)).^2 + (Y+(3*pitch*sin(4*pi/6))).^2 < (d/2).^2) = 1;
n((X+(3*pitch)).^2 + (Y+(4*pitch*sin(-5*pi/3))).^2 < (d/2).^2) = 1;

n((X-(pitch.*cos(5*pi/3))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(2*pi/3))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(3*pitch.*cos(2*pi/6))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(3*pitch.*cos(4*pi/6))).^2 + (Y+(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;

n((X-(pitch.*cos(5*pi/3))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(pitch.*cos(2*pi/3))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(3*pitch.*cos(2*pi/6))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;
n((X-(3*pitch.*cos(4*pi/6))).^2 + (Y-(5*pitch.*sin(5*pi/3))).^2 < (d/2).^2) = 1;

