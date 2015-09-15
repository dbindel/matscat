% plot_transmission(elt, Emin, Emax, npts)
%   Plot the transmission coefficient for scattering of a right-traveling
%   particle for energies Emin to Emax.  Also plot resonance poles with 
%   real part between Emin and Emax.
%
% Example:
%  elt = square_well([-5,-4,4,5], [0,5,0]);
%  plot_transmission(elt, 1, 10, 100);

function plot_transmission(elt, Emin, Emax, npts)

if nargin < 4, npts = 100; end

l = checked_resonances2(elt);
E = l.^2;

Es = linspace(Emin, Emax, npts);
ls = sqrt(Es);
[ts,rs] = compute_transmission(elt,ls);

subplot(2,1,1); 
plot(Es,abs(ts).^2); 
title('Transmission coefficient (top) and resonances (bottom)');
ylabel('Amplitude');
xlim([Es(1), Es(end)]);

subplot(2,1,2); 
E = E( find(imag(E) <= 0 & real(E) >= Emin & real(E) <= Emax) );
plot(real(E), imag(E), '*'); 
xlabel('Real energy');
ylabel('Imag energy');
xlim([Es(1), Es(end)]);
