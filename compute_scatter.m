% u = compute_scatter(elt,l)
%
% Compute the scattered wave in response to a plane wave of the form exp(i*l*x).

function u = compute_scatter(elt,l)

[K0,K1,K2] = form_operators(elt);
F = plane_forcing(elt,l);
u = (-l^2*K2 + 1i*l*K1 + K0)\F;
