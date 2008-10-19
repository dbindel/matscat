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
