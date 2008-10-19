function plot_resonances(elt, neigs)

if nargin < 2, neigs = 0; end

subplot(2,1,1)
plot_potential(elt);
title('Potential');

subplot(2,1,2)
%plot(checked_resonances(elt,neigs), '.', 'MarkerSize', 16);
plot(checked_resonances2(elt,neigs), '.', 'MarkerSize', 16);
title('Pole locations');

