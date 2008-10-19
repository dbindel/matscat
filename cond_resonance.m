% [kappa, Qpsi2, Qapsi2] = cond_resonance(elt,l,psi)
%
% Estimate a condition number for the resonance calculation.
%
% Input:
%   elt - mesh data structure
%   l   - the approximate eigenvalue
%   psi - the approximate eigenvector
%
% Output:
%   kappa -  the condition number with respect to max-norm perturbations in
%            the potential.  Equal to 
%             int |psi|^2 / ( 2*l * int psi^2 + i*(psi(a)^2 + psi(b)^2) )
%   Qpsi2  - the integral of psi^2 over (a,b)
%   Qapsi2 - the integral of |psi|^2 over (a,b)

function [kappa, Qpsi2, Qapsi2] = cond_resonance(elt,l,psi)

Qpsi2  = 0;  % Accumulate the integral of psi^2
Qapsi2 = 0;  % Accumulate the integral of |psi|^2

base = 1;
nelt = length(elt);
for j = 1:nelt
  
  % Element contribution (including matching conditions)
  order   = elt(j).order;
  [x,w]   = clencurt(order);
  w       = w * (elt(j).b - elt(j).a) / 2;
 
  % Assembly
  I             = base:base+order;
  Qpsi2         = Qpsi2  + w * psi(I).^2;
  Qapsi2        = Qapsi2 + w * abs(psi(I)).^2;

  base = base + order;
end

kappa = Qapsi2 / abs( 2*l*Qpsi2 + 1i*(psi(end)^2 + psi(1)^2) );
