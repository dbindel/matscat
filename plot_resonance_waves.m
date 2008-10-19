% plot_resonance_waves(elt, neigs, lrange)
%
% Plot the wave functions associated with bound states and resonances.
%
% Inputs:
%   elt    - Mesh data structure
%   neigs  - Number of eigenpairs to be computed (default: 0 -> compute all)
%   lrange - Pair of complex numbers [zmin, zmax] defining the lower
%            left and upper right corners of a box in the complex plane.
%            Eigenvalues and resonances in this box will be plotted
%            in order of descending imaginary part.  If no box is
%            provided, all computed eigenvalues and resonances are shown.
%   tol    - Absolute error tolerance for eigenvalues (default: 1e-6)

function [l,dl,V] = plot_resonance_waves(elt, neigs, lrange, tol)

if nargin < 2, neigs = 0;   end
if nargin < 3, lrange = []; end
if nargin < 4, tol = 1e-6;  end

[l,dl,V] = checked_resonances2(elt,neigs,tol);
if ~isempty(lrange)
  idx = find( real(l) >= real(lrange(1)) & imag(l) >= imag(lrange(1)) & ...
              real(l) <= real(lrange(2)) & imag(l) <= imag(lrange(2)) );
  l = l(idx);
  V = V(:,idx);
end

[junk,idx] = sort(-imag(l));
N = problem_size(elt);
for i = 1:length(idx)
  ii = idx(i);
  clf; plot_fields(elt, V(:,ii));
  title(sprintf('lambda = %g + %gi\n', real(l(ii)), imag(l(ii))));
  pause;
end

