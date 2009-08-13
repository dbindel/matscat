% [l] = checked_resonances2(elt,neigs,tol)
%
% Compute resonances with error estimate in order to check convergence.
% If the same answer occurs to within tol, accept the pole as converged.
% Inputs:
%   elt   - coarse mesh
%   neigs - number of poles desired? (default: 0 -> compute all)
%   tol   - absolute estimated error tolerance (default: 1e-6)
%           Use tol = 0 to return everything

function [l,dl,V] = checked_resonances2(elt,neigs,tol)

if nargin < 2, neigs = 0;    end
if nargin < 3, tol   = 1e-6; end

[l,V] = compute_resonances(elt, neigs);
if neigs == 0, neigs = length(l); end
N = problem_size(elt);
dl = zeros(neigs,1);
V = V(1:N,:);
for k = 1:neigs
  dl(k) = errest_resonance(elt, l(k), V(:,k));
end

if tol > 0  % Filter eigenvalues
  Igood = find(abs(dl) < tol);
  l = l(Igood);
  dl = dl(Igood);
  V = V(:,Igood);
end
