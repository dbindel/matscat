% elt = square_well(ab, V, h, order)
%
% Discretize a problem where the potential is given by a piecewise constant.
% Input:
%   ab    - breakpoints between regions where the potential is constant
%   V     - potential values on each interval ab(j) to ab(j+1)
%   h     - maximum element size (default: 1)
%   order - maximum element order (default: 20)

function elt = square_well(ab, V, h, order)

if nargin == 0
  ab = [-2, -1, 1, 2];
  V  = [0, -10, 0];
end

if nargin == 1
  V  = [0, ab, 0];
  ab = [-2, -1, 1, 2];
end

if nargin < 3, h = 1;      end
if nargin < 4, order = 20; end

elt = [];
nint = length(V);
i = 1;
for j = 1:nint
  L = ab(j+1)-ab(j);
  nelt = ceil(L/h);
  for k = 1:nelt
    elt(i).a = ab(j) + L*(k-1)/nelt;
    elt(i).b = ab(j) + L*k/nelt;
    elt(i).order = order;
    elt(i).Vtype = 'constant';
    elt(i).V = V(j);
    i = i + 1;
  end
end
