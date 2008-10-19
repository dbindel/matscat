%
% Compute resonances and error bounds for an Eckart potential
%
function [l,dl] = demo_eckart3

disp('Compute resonances for an Eckart potential barrier');
elt = func_well(@eckart_potential, -1:0.5:1, 20);
[l,dl] = checked_resonances2(elt, 30);


function Vx = eckart_potential(x)

Vx = (cosh(x)).^(-2);
