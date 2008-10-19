% plot_potential(elt)
%
% Plot the potential

function plot_potential(elt)

nelt = length(elt);
u = zeros(problem_size(elt),1);
base = 1;
for j = 1:nelt
  order = elt(j).order;
  x = cos(pi*(order:-1:0)/order)';
  xelt = elt(j).a * (1-x)/2 + elt(j).b * (1+x)/2;
  u(base:base+order) = eval_potential(elt(j), xelt);
  base = base + order + 1;
end
plot_fields(elt,u);
