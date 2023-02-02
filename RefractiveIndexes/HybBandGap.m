function n = HybBandGap(X, Y, FiberParams)

Clad = FiberParams.Dclad/2;
lambda = FiberParams.lambda;
r = FiberParams.dGe/2;
d = FiberParams.d;
pitch = FiberParams.pitch;
dn = FiberParams.dn;

n = ones([size(X)]);
n(X.^2+Y.^2<Clad.^2) = SilicaIndex(1064e-9);
for jk = 1:5
R = sqrt((X-(jk.*pitch)).^2 +(Y).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
R = sqrt((X+(jk.*pitch)).^2 +(Y).^2);
n(R<r) = SilicaIndex(lambda) + (dn)*(1 - (R(R < r)/r).^2);
end


% first ring
alpha = pi/3;
for jk = [1 2 4 5]
    n((X-(pitch.*cos(jk*alpha))).^2 + (Y-(pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
end


% 2nd ring
alpha = pi/6;
for jk = [1 2 3 4 5 7 8 9 10 11]
    if mod(jk,2)==0
        n((X-(2*pitch.*cos(jk*alpha))).^2 + (Y-(2*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
    else
        n((X-(sqrt(3)*pitch.*cos(jk*alpha))).^2 + (Y-(sqrt(3)*pitch.*sin(jk*alpha))).^2 < (d/2).^2) = 1;
    end
end

% 3rd ring
alpha = pi/3;

for jk = [1 2 4 5]
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

for jk = [1 2 4 5]
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

for jk = [1 2 4 5]
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

end