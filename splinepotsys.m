% function l = splinepotsys(VV1, VV2, R, xx, neigs)
%
% Compute resonances of a one-dimensional system of potentials with the
% form
%
%   -u'' + [V1(x), R(x); R(x), V2(x)]*u - lambda^2 u = 0
%
% where the potentials V1, V2, and R are given by splines with values
% value VV1(i), VV2(i), R(i) at each point xx(i).  The behavior at the
% end points is chosen according to the "not-a-knot" condition; if the
% end values are zero, the resulting (compactly supported) spline will
% be C^1 at the end points.
%
% Input:
%   VV1, VV2, R - Value of the potentials at the spline points
%   xx - Coordinates of the spline points
%   neigs - Max eigenvalues or resonances to be computed (default: 40)
%
% Output:
%   l - Vector of resolved resonance poles in the lambda (wave number) plane

function l = splinepotsys(VV1, VV2, R, xx, neigs);

if nargin < 5
  neigs = 40;
end

eltV1 = spline_well(xx,VV1);
eltV2 = spline_well(xx,VV2);
eltR  = spline_well(xx,R);

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
plot_splines(VV1,VV2,R,xx); 
title('Potential');

subplot(2,1,2);
plot(checked_resonances(elt), '.', 'MarkerSize', 25); 
title('Pole locations');


function plot_splines(VV1, VV2, R, xx)

xxx = linespace(min(xx), max(xx), 100);
pp1 = spline(xx, [0, VV1, 0], xxx);
pp2 = spline(xx, [0, VV2, 0], xxx);
ppR = spline(xx, [0, R,   0], xxx);
xh = ishold;
hold on
plot(xxx, pp1, 'r', 'LineWidth', 3);
plot(xxx, pp2, 'b', 'LineWidth', 3);
plot(xxx, ppR, 'k', 'LineWidth', 3, 'LineStyle', '-.');
if ~xh, hold off; end
