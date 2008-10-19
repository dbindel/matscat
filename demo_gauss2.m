%
% Compute resonances for a potential defined by a Gaussian function.
% This is pretty close to smooth.
%
function demo_gauss2

disp('Plot wave functions for bound states of a Gaussian well');
elt = func_well(@gauss_potential, -1:0.5:1, 30);
clf;
plot_resonance_waves(elt,20, [0, 10i]);

function Vx = gauss_potential(x)

Vx = -10*exp(-(3*x).^2);
