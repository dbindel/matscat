%
% Compute resonances for a potential defined by a Gaussian function.
% This is pretty close to smooth.
%
function demo_eckart2

disp('Compute resonances for an Eckart barrier');
for ll = 4:0.2:10
  fprintf('ll = %g\n', ll);
  elt = func_well(@eckart_potential, linspace(-ll,ll,ceil(2*ll)), 30);
  N   = problem_size(elt);
  [l,V] = compute_resonances(elt,10);
  [K0,K1,K2] = form_operators(elt,1);
  for j = 1:10
    [kappa, Qpsi2, Qapsi2] = cond_resonance(elt,l(j),V(1:N,j));
    KK = K0 - 1i*l(j)*K1 - l(j)^2*K2;
    R1  = KK*V(1:N,j);
    R2  = KK*V(N+1:N+N,j);
    fprintf('lambda = (%+g,%+g);\tkappa = %g;\trnorm = %g\n', ...
            real(l(j)), imag(l(j)), kappa, ...
            norm(R1) / norm(V(1:N,j)));
    %fprintf('  I[ psi ^2] = %g\n', abs(l(j)*Qpsi2));
    %fprintf('  I[|psi|^2] = %g\n', abs(Qapsi2));
    %fprintf('  End square = %g\n', abs(V(1,j)^2 + V(end,j)^2));
  end
  V(1:N,1) = V(1:N,1)/V(1,1);
  plot_fields(elt, V(1:N,1));
  axis([-10,10,-1,1]);
  pause;
end

function Vx = eckart_potential(x)

Vx = (cosh(x)).^(-2);
