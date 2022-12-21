function n = StepIndex(X, Y, fiberParams)

Dcore = fiberParams{1}; Dclad = fiberParams{2}; NA = fiberParams{3};
lambda = fiberParams{4};
alpha = SilicaLosses(lambda)*4.343;
kappa = (alpha.*log(10)/20)*(lambda)/(2*pi);
nClad = SilicaIndex(lambda);
nCore = sqrt(nClad^2 + NA^2);

n = ones(size(X));
n(X.^2+Y.^2<(Dclad/2).^2) = nClad;
n(X.^2+Y.^2<(Dcore/2).^2) = nCore;
n = n + 1i.*kappa;