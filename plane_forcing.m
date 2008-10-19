% F = plane_forcing(elt,l)
%
% Compute a forcing vector corresponding to the influence of an incident
% wave of the form exp(l*x)

function F = plane_forcing(elt,l)

z = l/1i;
F = zeros(problem_size(elt),1);
base = 1;
for j = 1:length(elt)
  order = elt(j).order;
  x = cos(pi*(order:-1:0)/order)';
  xelt = elt(j).a * (1-x)/2 + elt(j).b * (1+x)/2;
  F(base:base+order) = -eval_potential(elt(j), xelt) .* exp(z*xelt);
  F(base) = 0;
  F(base+order) = 0;
  base = base + order;
end
