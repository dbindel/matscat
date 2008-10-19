%
% Compute scattering states for a piecewise constant trapping potential
%

disp('Now we consider the scattering from a plane wave from a trapping');
disp('potential consisting of two short barriers.');
disp('Press any key to begin');
pause;

clf;
elt = square_well([-2,-1.1,-1,1,1.1,2], [0,10,0,10,0]);
subplot(2,1,1), plot_potential(elt); title('Potential');
subplot(2,1,2), 

kk = 1.7:0.05:2.3;
for j = 1:length(kk)
  title(sprintf('Scattered wave function (k = %g * pi)', kk(j)));
  u = compute_scatter(elt, kk(j)*pi);
  animate_wave(elt,u); 
end

