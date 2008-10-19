% elt = func_well(V, pts, order)
%
% Discretize a problem where the potential is given by a function handle
% Input:
%   V     - potential function (or function list)
%   pts   - element end points (default: -1:1)
%   order - maximum element order (default: 20)
%
% Additional arguments passed to func_well are passed on to the potential
% function V.

function elt = func_well(V, pts, order, varargin)

if nargin < 2, pts = -1:1;  end
if nargin < 3, order = 20;  end

x = cos(pi*(order:-1:0)/order)';
for i = 1:length(pts)-1
  elt(i).a = pts(i);
  elt(i).b = pts(i+1);
  elt(i).order = order;
  if iscell(V)
    elt(i).V = V{i};
  else
    elt(i).V = V;
  end
  elt(i).Vtype = 'function';
  elt(i).args = varargin;
end
