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
