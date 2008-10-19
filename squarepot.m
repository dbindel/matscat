% function l = squarepot(VV, xx, neigs)
%
% Compute resonances of a one-dimensional piecewise-constant potential.
% The value of the potential will be VV(i) on the interval ( xx(i), xx(i+1) ).
%
% Input:
%   VV - Value of the potential on each interval
%   xx - Coordinates of the interval endpoints
%   neigs - Max eigenvalues or resonances to be computed (default: 40)
%
% Output:
%   l - Vector of resolved resonance poles in the lambda (wave number) plane

function l = squarepot(VV,xx, neigs);

if nargin ==0
  VV = [0, 1, 0];
  xx = [-2,-1,1,2];
end

if nargin < 3
  neigs = 40;
end

elt = square_well(xx,VV);
l = checked_resonances(elt,neigs);

clf;
subplot(2,1,1) ;
plot_potential1(VV,xx); 
title('Potential');

subplot(2,1,2);
plot(checked_resonances(elt), '.', 'MarkerSize', 25); 
title('Pole locations');

