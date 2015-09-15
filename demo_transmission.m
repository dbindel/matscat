elt = square_well([-5,-4,4,5], [0,5,0]);

clf;
disp('First, show the resonances for a square barrier');
plot_resonance(elt);
pause;

clf;
disp('Now show the transmission coefficient and resonances in the E plane');
plot_transmission(elt, 1, 10, 500);
pause;
