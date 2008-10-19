function [Dmf2,f2] = chebdifft(f,M,Nout)

% The function [Df2,f2] = chebdifft(f,M,Nout) computes the M'th
% approximate Chebyshev derivatives of the data vector y and interpolates
% both the derivative and the original data onto a different Chebyshev grid.
% A Fast Fourier Transform is used compute the Chebyshev cofficients
% of the data vector. A recursion formula is used to compute the
% Chebyshev coefficients for each derivative. A FFT is then used again
% to compute the derivatives in physical space.
% 
%  Input:
%  f:        Vector containing function values at the Chebyshev points
%            x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
%  M:        Derivative required (positive integer)
%  Nout:     Number of points in the output mesh (may be different from input)
%
%  Output:
%  Df2:      Vector containing approximate M'th derivative at
%            y(k) = cos((k-1)*pi/(Nout-1)), k = 1...Nout.
%  f2:       Vector containing interpolated values at
%            y(k) = cos((k-1)*pi/(Nout-1)), k = 1...Nout.
%            
%  J.A.C. Weideman, S.C. Reddy 2000.  Modified dbindel, 2006.


f=f(:);                                      % Make sure f is a vector     
N=length(f);      
a0=fft([f; flipud(f(2:N-1))]);               % Extend and compute fft
a0=a0(1:N).*[0.5; ones(N-2,1); 0.5]/(N-1);   % a0 contains Chebyshev 
                                             % coefficients of f

a=[a0 zeros(N,M)];                           % Recursion formula
for ell=1:M                                  % for computing coefficients
  a(N-ell,ell+1)=2*(N-ell)*a(N-ell+1,ell);   % of ell'th derivative 
  for k=N-ell-2:-1:1
    a(k+1,ell+1)=a(k+3,ell+1)+2*(k+1)*a(k+2,ell);
  end;
  a(1,ell+1)=a(2,ell)+a(3,ell+1)/2;
end;

if Nout < N
  a = a(1:Nout,:);
else
  a = [a; zeros(Nout-N,M+1)];
end

dback=[2*a(1,M+1); a(2:Nout-1,M+1); 2*a(Nout,M+1); flipud(a(2:Nout-1,M+1))];
back =[2*a(1,  1); a(2:Nout-1,  1); 2*a(Nout,  1); flipud(a(2:Nout-1,  1))];
Dmf2=0.5*fft(dback);                           % Transform back to
Dmf2=Dmf2(1:Nout);                             % physical space
f2=0.5*fft(back);
f2=f2(1:Nout);

if max(abs(imag(f))) == 0    % Real data in, real data out
  Dmf2 = real(Dmf2); 
  f2   = real(f2);
end
