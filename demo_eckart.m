%
% Compute resonances for a potential defined by a Gaussian function.
% This is pretty close to smooth.
%
function demo_eckart

disp('Compute resonances for an Eckart barrier');
%for ll = 4:0.2:10
  ll = 6;
  elt = func_well(@eckart_potential, linspace(-ll,ll,ceil(2*ll)), 20);
  clf; 
  plot_resonance(elt,20); 
  subplot(2,1,1); hold on; axis([-10,10,0,1.2]);
  subplot(2,1,2); hold on; axis([-2,2,-2,0]);
  pause(0.2);
%end

function Vx = eckart_potential(x)

Vx = (cosh(x)).^(-2);
