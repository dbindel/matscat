% ts = compute_transmission(elt,l)
%
% Compute the transmission in response to a plane wave of the form exp(i*l*x).

function [ts,rs] = compute_scatter(elt,ls)

ts = 0*ls;
rs = 0*ls;
xx = get_xgrid(elt);
for j = 1:length(ls)
  l = ls(j);
  [K0,K1,K2] = form_operators(elt);
  F = plane_forcing(elt,l);
  u = (-l^2*K2 + 1i*l*K1 + K0)\F;
  ts(j) = 1 + u(end) / exp(-1i*l*xx(end));
  rs(j) = u(1)   / exp( 1i*l*xx(1));
end
