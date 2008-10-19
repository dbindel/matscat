%
% Compute scattering states and resonances for a one-dimensional square well
% potential ("particle in a box").
%

% -- Example 1: Show dynamics of resonances changing to bound states --

disp('In the first example, we show how the eigenvalues and resonances');
disp('change as the depth of a square potential well changes.');
disp('Press any key to begin');
pause;

clf;
V0s = linspace(6,10, 20);
for jj = 1:length(V0s)
  elt = square_well(-V0s(jj));   % Square well on [-1,1] with depth V0
  plot_resonance(elt,20);        % Plot potential and resonance
  pause(0.1);
end

% -- Example 2: Show scattering --

disp(' ');
disp('Now we consider the scattering from a plane wave of the form');
disp('exp(-ikx), where k = pi.  The real part of the scattered wave');
disp('is in red; the imaginary part is in blue.');
disp('Press any key to begin');
pause;

clf;
elt = square_well(-10);
u = compute_scatter(elt, -pi);
subplot(2,1,1), plot_potential(elt); title('Potential');
subplot(2,1,2), plot_fields(elt,u);  title('Scattered wave function');


% -- Example 3: Animation --

disp(' ');
disp('In addition to just displaying the potential, we can also animate');
disp('the wave -- that is, we show Re(exp(i*t)*u_s), where u_s is the');
disp('scattered wave function.  We will show an animation for');
disp('a range of wave numbers for the incident wave');
disp('Press any key to begin');
pause;

ks = 0.8:0.2:1.2;
for jj = 1:length(ks)
  title(sprintf('Scattering from exp(i * %g * pi * x)', ks(jj)));
  u_s = compute_scatter(elt, ks(jj)*pi);
  animate_wave(elt, u_s);
end
