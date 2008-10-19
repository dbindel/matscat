% elt = spline_well(xx, Vx, h, order)
%
% Discretize a problem where the potential is given by a spline.
% Input:
%   xx    - interpolation points
%   Vx    - spline values at the points
%   h     - maximum element size (default: 1)
%   order - maximum element order (default: 20)

function elt = spline_well(xx, Vx, h, order)

if nargin < 3, h = 1;      end
if nargin < 4, order = 20; end

elt = [];
L = xx(end)-xx(1);
pp = spline(xx,Vx);
nelt = ceil(L/h);
for i = 1:nelt
  elt(i).a = xx(1) + L*(i-1)/nelt;
  elt(i).b = xx(1) + L*i/nelt;
  elt(i).order = order;
  elt(i).Vtype = 'spline';
  elt(i).V = pp;
end
