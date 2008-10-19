% animate_wave(elt,u,N)
%
% Animate one period of the wave solution.
% Inputs:
%   elt - mesh description
%   u   - complex solution vector (e.g. from compute_scatter)
%   N   - number of cycles in the animation 
%
function animate_wave(elt,u,N)

if nargin < 4, N = 24; end
umax = max(abs(u));
for j = 1:N
  plot_fields(elt, real(u*exp(j*2i*pi/N)));
  axis([elt(1).a, elt(end).b, -umax, umax]);
  pause(0.2);
end
