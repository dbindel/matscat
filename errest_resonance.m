% [dl] = errest_resonance(elt, l, psi)
%
% Estimate the sensitivity of the resonance calculation.
% Input:
%   elt - Mesh data structure
%   l   - Computed eigenvalue
%   psi - Computed wave function
%
% Output:
%   dl  - Approximate error computed from linearized perturbation theory

% FIXME: I ought to make the chebdifft2 indexing convention equivalent
% (and the clencurt convention) equivalent to the convention for cheb

function [dl] = errest_resonance(elt, l, psi)

Qpsi2  = 0;  % Accumulate the integral of psi^2
QpsiR  = 0;  % Accumulate integral of psi*R, R = (H-l^2) psi

base = 1;
nelt = length(elt);
for j = 1:nelt
  
  % Element contribution (including matching conditions)
  order   = elt(j).order;
  [xf,wf] = clencurt(5*order); xf = xf(end:-1:1);
  xf      = elt(j).a * (1-xf)/2 + elt(j).b * (1+xf)/2;
  wf      = wf * (elt(j).b - elt(j).a) / 2;
  Vf      = eval_potential(elt(j), xf);
 
  % Assembly
  I             = base:base+order;
  [dpsif, psif] = chebdifft2(psi(I), 2, 5*order+1);
  dpsif         = dpsif * 4 / (elt(j).b - elt(j).a)^2;
  Rf            = -dpsif + Vf.*psif - l^2*psif;

  Qpsi2         = Qpsi2 + wf * psif.^2;
  QpsiR         = QpsiR + wf * (psif .* Rf);

  base = base + order;
end

dl = QpsiR / (2*l*Qpsi2 + 1i*(psi(1)^2 + psi(end)^2));
