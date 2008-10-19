%
% Compute scattering states and resonances for a potential defined by
% spline function.  For this problem, the potential is C^1.
%

disp('Compute scattering from a function defined by a spline');

clf;
V0s = linspace(0, 12, 24);
for jj = 1:length(V0s)
  elt = spline_well([-2,-1,1,2], [0, 0,-V0s(jj),-V0s(jj),0, 0]);
  subplot(2,1,1), plot_potential(elt);                  axis([-2,2, -20,0]);
  subplot(2,1,2), plot(checked_resonances2(elt), '.');  axis([-5,5, -3,3]);
  pause(0.2);
end
