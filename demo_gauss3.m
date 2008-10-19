%
% Compute resonances for a potential defined by a Gaussian function.
% This is pretty close to smooth.
%
function [l,dl] = demo_gauss3

disp('Compute resonances and error bounds for a Gaussian potential well');
elt = func_well(@gauss_potential, -1:0.5:1, 20);
[l,dl] = checked_resonances2(elt, 30, 0);


function Vx = gauss_potential(x)

Vx = -10*exp(-(3*x).^2);
