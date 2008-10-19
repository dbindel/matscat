%
% Compute resonances for a potential defined by a Gaussian function.
% This is pretty close to smooth.
%
function demo_Vwell

disp('Compute resonances for barrier trapping');
for a = 20:2:80
  elt = func_well({@V_potential, @well_potential, ...
                   @well_potential, @V_potential}, -2:1:2, 20, a);
  plot_resonance(elt, 40);
  shg;
  pause(0.5);
end


function Vx = V_potential(x,a)

x = abs(x);
Vx = a*max(1 - 2*abs(x-1.5), 0);


function Vx = well_potential(x,a)

Vx = 0*x - a;
