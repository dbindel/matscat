% plot_potential1(VV1,VV2,R,xx)
%
% Plot the potential with continuous lines

function plot_potential1s(VV1,VV2,R,xx)

VV1 = [0,VV1,0];
VV2 = [0,VV2,0];
R   = [0,R,  0];

xx = [2*xx(1)-xx(2), xx, 2*xx(end)-xx(end-1)];
yy = zeros(2*length(xx),1);
yy(1:2:end-1) = xx;
yy(2:2:end)   = xx;
WW1 = zeros(2*length(xx),1);
WW1(2:2:end-2) = VV1;
WW1(3:2:end-1) = VV1;
WW2 = zeros(2*length(xx),1);
WW2(2:2:end-2) = VV2;
WW2(3:2:end-1) = VV2;
WW3 = zeros(2*length(xx),1);
WW3(2:2:end-2) = R;
WW3(3:2:end-1) = R;

xh = ishold;
hold on
plot(yy,WW1,'r','LineWidth',3);
plot(yy,WW2,'b','LineWidth',3);
plot(yy,WW3,'k','LineWidth',3, 'LineStyle', '-.');
if ~xh, hold off; end

VV = [VV1,VV2,R];
ylim([min(VV)-(max(VV)-min(VV))/3, max(VV)+(max(VV)-min(VV))/3])
