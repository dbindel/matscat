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
% [l] = checked_resonances(elt,neigs,tol)
%
% Compute resonances with two densities in order to check convergence.
% If the same answer occurs to within tol, accept the pole as converged.
% Inputs:
%   elt   - coarse mesh
%   neigs - number of poles desired? (default: 0 -> compute all)
%   tol   - match tolerance (default: 1e-6)

function l = checked_resonances(elt,neigs,tol)

if nargin < 2, neigs = 0;  end
if nargin < 3, tol = 1e-6; end

elt2 = elt;
for j = 1:length(elt)
  elt2(j).order = ceil(elt2(j).order*1.5);
end
l1 = compute_resonances(elt,  neigs);
l2 = compute_resonances(elt2, neigs);

l = compare_eigs(l2,l1,tol);
% [l,V] = compute_resonances(elt, neigs)
%
% Compute resonances via a generalized linear eigenvalue problem.
% Inputs:
%   elt   - Problem description
%   neigs - Number of poles desired? (default: 0 --> compute all)

function [l,V] = compute_resonances(K0,K1,K2,neigs)

if nargin == 1, neigs = 0;  end
if nargin == 2, neigs = K1; end
if nargin == 3, neigs = 0;  end

if nargin == 1 | nargin == 2
  elt = K0;
  [N,nnz] = problem_size(elt);
  [K0,K1,K2] = form_operators(elt, ...
     (neigs ~= 0) & (nnz < 0.2 * N^2) & (N > 100));
end

N = length(K0);
if neigs == 0 | issparse(K0) ~= 0
  Z = zeros(N);
  I = eye(N);
else
  Z = spalloc(N,N,0);
  I = speye(N);
end
A = [K0, Z; Z, I];
B = [-K1, -K2; I, Z];

if neigs == 0
  if nargout == 1
    l = eig(A,B)*1i;
  elseif nargout == 2
    [V,D] = eig(A,B);
    l = diag(D)*1i;
  end
elseif issparse(K0)
  [L,U,P,Q] = lu(A);
  opts.isreal = 1;
  opts.disp = 0;
  if nargout == 1
    l = eigs(@cr_lufun1, 2*N, neigs, 0, opts, L, U, P, Q, B)*1i;
  elseif nargout == 2
    [V,D] = eigs(@cr_lufun1, 2*N, neigs, 0, opts, L, U, P, Q, B);
    l = diag(D)*1i;
  end
else
  [L,U] = lu(A);
  opts.isreal = 1;
  opts.disp = 0;
  if nargout == 1
    l = eigs(@cr_lufun2, 2*N, neigs, 0, opts, L, U, B)*1i;
  elseif nargout == 2
    [V,D] = eigs(@cr_lufun2, 2*N, neigs, 0, opts, L, U, B);
    l = diag(D)*1i;
  end
end

function y = cr_lufun1(x,L,U,P,Q,B)
y = Q*(U\(L\(P*(B*x))));

function y = cr_lufun2(x,L,U,B)
y = U\(L\(B*x));
% [ll,diff] = compare_eigs(ll1, ll2, tol)
%
% Compare the results of two eigenvalue vectors (ll1 and ll2).
% Inputs:
%   ll1, ll2 - eigenvalue lists
%   tol      - matching tolerance
% Outputs:
%   ll   - all elements of ll1 within tol of some element of ll2
%   diff - diff(i) is min(abs(ll1(i) - ll2));

function [ll,diff] = compare_eigs(ll1, ll2, tol)

if size(ll1,1) == 1, ll1 = ll1.'; end
if size(ll2,1) ~= 1, ll2 = ll2.'; end
e1 = 0*ll1 + 1;
e2 = 0*ll2 + 1;
diff = min(abs(ll1*e2 - e1*ll2), [], 2);
ll   = ll1(find(diff < tol));
% plot_potential1(VV,xx)
%
% Plot the potential with continuous lines

function plot_potential1(VV,xx)

VV = [0,VV,0];
xx = [2*xx(1)-xx(2), xx, 2*xx(end)-xx(end-1)];
yy = zeros(2*length(xx),1);
WW = zeros(2*length(xx),1);
yy(1:2:end-1) = xx;
yy(2:2:end)   = xx;
WW(2:2:end-2) = VV;
WW(3:2:end-1) = VV;

plot(yy,WW,'r','LineWidth',3)
ylim([min(VV)-(max(VV)-min(VV))/3, max(VV)+(max(VV)-min(VV))/3])
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
% [K0,K1,K2] = form_operators(elt,is_sparse)
%
% Form pseudospectral discretization of
%   ((-D^2 + V) + z^2) u = 0  on  (a,b)
%   (D-z) u = 0               at  a
%   (D+z) u = 0               at  b
% The matrices K0, K1, K2 are such that
%   K(l) = K0 + z*K1 + z^2*K2
% is the discretized operator.

function [K0,K1,K2] = form_operators(elt,is_sparse)

if nargin < 2, is_sparse = 0; end

if ~isfield(elt, 'V')
  [K0,K1,K2] = form_operators_sys(elt, is_sparse);
  return
end

% Useful constants and memory allocations
nelt  = length(elt);
if is_sparse
  [N,nz] = problem_size(elt);
  K0 = spalloc(N,N,nz);
  K1 = spalloc(N,N,2);
  K2 = spalloc(N,N,N);
else
  N  = problem_size(elt);
  K0 = zeros(N);
  K1 = zeros(N);
  K2 = zeros(N);
end

% Operator on subdomains
base = 1;
for j = 1:nelt
  
  % Element contribution (including matching conditions)
  order   = elt(j).order;
  [D,x]   = cheb(order);
  lelt(j) = elt(j).b - elt(j).a;
  xelt    = elt(j).a * (1-x)/2 + elt(j).b * (1+x)/2;
  K0elt   = -4*D^2/lelt(j)^2 + diag(eval_potential(elt(j), xelt));
  K0elt(1,:)   = -2*D(1,:)/lelt(j);
  K0elt(end,:) =  2*D(end,:)/lelt(j);
  
  % Assembly
  I             = base:base+order;
  Iint          = base+1:base+order-1;
  K0(I,I)       = K0(I,I) + K0elt;
  K2(Iint,Iint) = eye(order-1);
  
  base = base + order;
end

% Radiation boundary conditions
K1(1,1) = 1;
K1(N,N) = 1;


%
% Compute D = differentiation matrix, x = Chebyshev grid.  Adapted
% from the routine in Trefethen's "Spectral Elements in MATLAB"
%
function [D,x] = cheb(N)

x = cos(pi*(0:N)/N)'; 
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';                  
D  = (c*(1./c)')./(dX+(eye(N+1)));      % off-diagonal entries
D  = D - diag(sum(D'));                 % diagonal entries
x = x(N+1:-1:1);
D = D(N+1:-1:1,N+1:-1:1);
% N = problem_size(elt)
%
% Get total number of points in the mesh

function [N,nz] = problem_size(elt)

N = 1;
for k = 1:length(elt)
  N = N + elt(k).order;
end

if nargout == 2
  nz = 0;
  for k = 1:length(elt)
    nz = nz + (elt(k).order+1)^2;
  end
end
% [Vx] = eval_potential(eltj, x)
%
% Evaluate the potential associated with the mesh element eltj at x.

function [Vx] = eval_potential(eltj, x)

if strcmp(eltj.Vtype, 'function')
  Vx = feval(eltj.V, x, eltj.args{:});
elseif strcmp(eltj.Vtype, 'spline')
  Vx = ppval(eltj.V, x);
elseif strcmp(eltj.Vtype, 'constant')
  Vx = eltj.V + 0*x;
else
  error('Unknown potential type');
end
% [K0,K1,K2] = form_operators_sys(elt, is_sparse)
%
% Form pseudospectral discretization of
%   ((-D^2 + V1) + z^2) u + R v = 0 on (a,b)
%   ((-D^2 + V2) + z^2) v + R u = 0 on (a,b)
%   (D-z) u = (D-z) v = 0           at  a
%   (D+z) u = (D+z) v = 0           at  b
% The matrices K0, K1, K2 are such that
%   K(l) = K0 + z*K1 + z^2*K2
% is the discretized operator.

% FIXME: The setup for this ought to be cleaner...

function [K0,K1,K2] = form_operators_sys(elt, is_sparse)

if nargin < 3, is_sparse = 0; end

[K0a, K1a, K2a] = form_operators(extract_eltV1(elt), is_sparse);
[K0b, K1b, K2b] = form_operators(extract_eltV2(elt), is_sparse);

% Useful constants and memory allocations
eltR = extract_eltR(elt);
nelt = length(eltR);
if is_sparse
  [N,nz] = problem_size(eltR);
  K0R = spalloc(N,N,nz);
  Z   = spalloc(N,N,nz);
else
  N   = problem_size(elt);
  K0R = zeros(N);
  Z   = zeros(N);
end

% Operator on subdomains
base = 1;
for j = 1:nelt
  
  % Element contribution (including matching conditions)
  order = elt(j).order;
  x     = cos(pi*(order:-1:0)/order)'; 

  lelt(j) = elt(j).b - elt(j).a;
  xelt    = elt(j).a * (1-x)/2 + elt(j).b * (1+x)/2;
  K0elt   = diag(eval_potential(eltR(j), xelt));
  K0elt(1,:)   = 0;
  K0elt(end,:) = 0;

  % Assembly
  I             = base:base+order;
  Iint          = base+1:base+order-1;
  K0R(I,I)      = K0R(I,I) + K0elt;
  
  base = base + order;
end

K0 = [K0a, K0R; K0R, K0b];
K1 = [K1a, Z;   Z,   K1b];
K2 = [K2a, Z;   Z,   K2b];

% plot_potential1(VV1,VV2,R,xx)
%
% Plot the potential with continuous lines

function plot_potential1s(VV1,VV2,R,xx)

VV1 = [0,VV1,0];
VV2 = [0,VV2,0];
R   = [0,R,  0];

xx = [2*xx(1)-xx(2), xx, 2*xx(end)-xx(end-1)];
yy = zeros(2*length(xx),1);
yy(1:2:end-1) = xx;
yy(2:2:end)   = xx;
WW1 = zeros(2*length(xx),1);
WW1(2:2:end-2) = VV1;
WW1(3:2:end-1) = VV1;
WW2 = zeros(2*length(xx),1);
WW2(2:2:end-2) = VV2;
WW2(3:2:end-1) = VV2;
WW3 = zeros(2*length(xx),1);
WW3(2:2:end-2) = R;
WW3(3:2:end-1) = R;

xh = ishold;
hold on
plot(yy,WW1,'r','LineWidth',3);
plot(yy,WW2,'b','LineWidth',3);
plot(yy,WW3,'k','LineWidth',3, 'LineStyle', '-.');
if ~xh, hold off; end

VV = [VV1,VV2,R];
ylim([min(VV)-(max(VV)-min(VV))/3, max(VV)+(max(VV)-min(VV))/3])
function [elt] = extract_eltV1(elt)

for k = 1:length(elt)
  elt(k).V     = elt(k).V1;
  elt(k).Vtype = elt(k).V1type;
end
function [elt] = extract_eltV2(elt)

for k = 1:length(elt)
  elt(k).V     = elt(k).V2;
  elt(k).Vtype = elt(k).V2type;
end
function [elt] = extract_eltR(elt)

for k = 1:length(elt)
  elt(k).V     = elt(k).R;
  elt(k).Vtype = elt(k).Rtype;
end
