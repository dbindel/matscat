% plot_potential1(VV,xx)
%
% Plot the potential with continuous lines

function plot_potential1(VV,xx)

VV = [0,VV,0];
xx = [2*xx(1)-xx(2), xx, 2*xx(end)-xx(end-1)];
yy = zeros(2*length(xx),1);
WW = zeros(2*length(xx),1);
yy(1:2:end-1) = xx;
yy(2:2:end)   = xx;
WW(2:2:end-2) = VV;
WW(3:2:end-1) = VV;

plot(yy,WW,'r','LineWidth',3)
ylim([min(VV)-(max(VV)-min(VV))/3, max(VV)+(max(VV)-min(VV))/3])
