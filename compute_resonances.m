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

if nargin == 1 || nargin == 2
  elt = K0;
  [N,nnz] = problem_size(elt);
  [K0,K1,K2] = form_operators(elt, ...
     (neigs ~= 0) & (nnz < 0.2 * N^2) & (N > 100));
end

N = length(K0);
if neigs == 0 || issparse(K0) ~= 0
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

  cr_lufun1 = @(L, U, P, Q, B, x)( @(x)( Q*(U\(L\(P*(B*x)))) ) );

  if nargout == 1
    l = eigs(cr_lufun1(L,U,P,Q,B), 2*N, neigs, 0, opts)*1i;
  elseif nargout == 2
    [V,D] = eigs(cr_lufun1(L, U, P, Q, B), 2*N, neigs, 0, opts);
    l = diag(D)*1i;
  end
else
  [L,U] = lu(A);
  opts.isreal = 1;
  opts.disp = 0;

  cr_lufun2 = @(L, U, B, x)( @(x)( U\(L\(B*x)) ) );

  if nargout == 1
    l = eigs(cr_lufun2(L,U,B), 2*N, neigs, 0, opts)*1i;
  elseif nargout == 2
    [V,D] = eigs(cr_lufun2(L,U,B), 2*N, neigs, 0, opts);
    l = diag(D)*1i;
  end
end
