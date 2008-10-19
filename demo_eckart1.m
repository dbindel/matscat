%
% Compute resonances for an Eckart barrier
%
function demo_eckart1

disp('Compute resonances for an Eckart barrier with varying cutoff');
for a = 0.1:0.1:3
  elt = func_well(@eckart_potential, -8:1:8, 30, a);
  clf;
  hold on
  plot(compute_resonances(elt,40), 'g*');
  plot(checked_resonances2(elt,40), 'b.');
%  plot(checked_resonances(elt,40), 'b.');
  mm = [];
  nu = roots([1,1,a^2]);
  for k=1:3
    mm = [mm, 1i*nu - 1i*(k-1)];
  end
  plot(mm,'or')
  pause(0.5);
  hold off
end

function Vx = eckart_potential(x,a)

Vx = a^2*sech(x).^2;
