% function l = squarepotsys(VV1, VV2, R, xx, neigs)
%
% Compute resonances of a one-dimensional system of potentials with the
% form
%
%   -u'' + [V1(x), R(x); R(x), V2(x)]*u - lambda^2 u = 0
%
% where the potentials V1, V2, and R are piecewise constant with values
% VV1(i), VV2(i), R(i) on (xx(i), xx(i+1)).
%
% Input:
%   VV1, VV2, R - Value of the potentials at the spline points
%   xx - Coordinates of the spline points
%   neigs - Max eigenvalues or resonances to be computed (default: 40)
%
% Output:
%   l - Vector of resolved resonance poles in the lambda (wave number) plane

function l = squarepot(VV1, VV2, R, xx, neigs);

% Compute resonances for a one-dimensional potential consisting of 
% number of steps
% The data is a (N-1) vector VV of values at of the potential
% between the values of an N vector  xx

if nargin < 5
  neigs = 40;
end

eltV1 = square_well(xx,VV1);
eltV2 = square_well(xx,VV2);
eltR  = square_well(xx,R);

elt = eltV1;
elt = rmfield(elt, 'V');
elt = rmfield(elt, 'Vtype');
for k = 1:length(eltV1)
  elt(k).V1 = eltV1(k).V; elt(k).V1type = eltV1(k).Vtype;
  elt(k).V2 = eltV2(k).V; elt(k).V2type = eltV2(k).Vtype;
  elt(k).R  = eltR(k).V;  elt(k).Rtype  = eltR(k).Vtype;
end

l = checked_resonances(elt,neigs);

clf;
subplot(2,1,1) ;
plot_potential1s(VV1,VV2,R,xx); 
title('Potential');

subplot(2,1,2);
plot(checked_resonances(elt), '.', 'MarkerSize', 25); 
title('Pole locations');

