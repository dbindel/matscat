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
