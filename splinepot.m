% function l = splinepot(VV, xx, neigs)
%
% Compute resonances of a one-dimensional potential obtained by taking
% a spline with value VV(i) at each point xx(i).  The behavior at the
% end points is chosen according to the "not-a-knot" condition; if the
% end values are zero, the resulting (compactly supported) spline will
% be C^1 at the end points.
%
% Input:
%   VV - Value of the potential at the spline points
%   xx - Coordinates of the spline points
%   neigs - Max eigenvalues or resonances to be computed (default: 40)
%
% Output:
%   l - Vector of resolved resonance poles in the lambda (wave number) plane

function l = splinepot(VV,xx, neigs);

% Default spline
if nargin ==0
  VV = [0,0.5,1,0.5, 0];
  xx = [-1,-0.5,0,0.5,1];
end

% By default compute the 40 closest to zero
if nargin < 3
  neigs = 40;
end

% Add zero knots at the end if length(VV) < length(xx)
if length(xx) > length(VV)
  VV = [0,VV,0];
end

elt = spline_well(xx,[0,VV,0]);
l = checked_resonances(elt,neigs);

clf;
subplot(2,1,1), 
plot_spline(xx,VV); 
title('Potential');

subplot(2,1,2), 
plot(checked_resonances(elt), '.', 'MarkerSize', 25);
title('Pole locations');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_spline(xx,VV)

[nn,NN]=size(xx);
xxx = linspace(max(xx),min(xx),100);
pp = spline(xx,[0,VV,0],xxx);
plot(xxx,pp,'r','LineWidth',3);
ylim([min(VV)-(max(VV)-min(VV))/3,max(VV)+(max(VV)-min(VV))/3]);
