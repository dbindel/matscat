%
% Compute resonances for a potential defined by a Gaussian function.
% This is pretty close to smooth.
%
function demo_gauss

disp('Compute resonances for a Gaussian potential well');
elt = func_well(@gauss_potential, -1:0.5:1, 30);
clf;
plot_resonance(elt,20);


function Vx = gauss_potential(x)

Vx = -10*exp(-(3*x).^2);
