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

