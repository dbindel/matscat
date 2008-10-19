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
