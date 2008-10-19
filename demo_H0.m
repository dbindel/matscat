%
% Compute scattering states and resonances for a potential defined by
% spline function.  For this problem, the potential is C^1.
%

disp('Compute resonance poles for a zero potential.  Show all the eigenvalues');
disp('of the matrix problem, with a special mark for the eigenvalues that');
disp('stay put under refinement.  The zero eigenvalue is (appropriately) the');
disp('only one marked as correct');

elt = spline_well([0,1], [0,0], 0.5, 20);
l1 = compute_resonances(elt);
l2 = checked_resonances2(elt);

clf;
plot(real(l1), imag(l1), 'b.'); hold on
plot(real(l2), imag(l2), 'ro');
legend('All eigenvalues', 'Checked eigenvalues');
