%
% Compute resonances for a potential defined by an Eckart potential *
% a continuous cutoff function defined by a spline.
%
function demo_eckart2

ll = 20;
elt = func_well(@eckart_potential, linspace(-ll,ll,ceil(4*ll)), 40, ll);
clf; 
plot(compute_resonances(elt,40), 'b.');
%plot(checked_resonances2(elt,40), 'b.');


function Vx = eckart_potential(x,ll)

sp1 = spline([-ll,-ll+3], [0,0, 1,0]);
sp2 = spline([ ll-3, ll], [1,0, 0,0]);
I1 = find(x < -ll+3);
I2 = find(x >  ll-3);
wt = 1+0*x;
if ~isempty(I1), wt(I1) = ppval(sp1, x(I1)); end
if ~isempty(I2), wt(I2) = ppval(sp2, x(I2)); end
Vx = (cosh(0.25*x)).^(-2) .* wt;
