% [xx] = get_xgrid(elt);
%
% Get the grid points associated with the spectral element mesh.

function x = get_xgrid(elt)

N = problem_size(elt);
x = zeros(N,1);

base = 1;
for j = 1:length(elt)
  order = elt(j).order;
  xc   = cos(pi*(order:-1:0)/order)';
  xelt = elt(j).a * (1-xc)/2 + elt(j).b * (1+xc)/2;
  x(base:base+order) = xelt;
  base = base + order;
end
