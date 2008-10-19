% plot_fields(elt,u)
%
% Plot some function defined on the spectral element mesh

function plot_fields(elt,u)

rflag = isreal(u);
hflag = ishold;
N = problem_size(elt);

hold on;
cla;
base = 1;
warning('OFF', 'MATLAB:polyfit:RepeatedPoints');
for j = 1:length(elt)
  order = elt(j).order;
  x = cos(pi*(order:-1:0)/order)';

  % Plot points  
  xelt = elt(j).a * (1-x)/2 + elt(j).b * (1+x)/2;
  uelt = u(base:base+order);
  plot(xelt, real(uelt), 'r.', 'markersize', 16);
  if ~rflag, plot(xelt, imag(uelt), 'b.', 'markersize', 16); end
  
  % Plot interpolant
  xxelt = linspace(elt(j).a, elt(j).b);
  [pc,s,mu] = polyfit(xelt,uelt,order);
  uuelt = polyval(pc, (xxelt-mu(1))/mu(2));
  hr = line(xxelt, real(uuelt));
  set(hr, 'Color', 'r');
  if ~rflag
    hi = line(xxelt, imag(uuelt));
    set(hi, 'Color', 'b'); 
  end

  if length(u) > N
    base = base + order + 1;
  else
    base = base + order;
  end
end
if ~ishold, hold off; end
warning('ON', 'MATLAB:polyfit:RepeatedPoints');
