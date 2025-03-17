load('results_fig5-6.mat')
gg=grid_fine;
R=linspace(0.1,0.9,grid_fine);
for i=1:gg
    diam2(i)=(32*f(i)-36*f(i)^2+9*f(i)^3)/(9*f(i)^2)-(4/9)*((-(-64*f(i)^2+144*f(i)^3-108*f(i)^4+27*f(i)^5)/(f(i)^4)))^(1/2);
    diam3(i)=(1/36)*((f(i)*(-5+3*f(i))*(-4+3*f(i)))/(-1+f(i))^2-2*(-(f(i)^2*(-4+3*f(i))^3)/(-1+f(i))^4)^(1/2));
    diam4(i)=(f(i)*(-4+3*f(i))*(-40+27*f(i)))/(32-27*f(i))^2-12*(-(f(i)^2*(-4+3*f(i))^3)/(-32+27*f(i))^4)^(1/2);
    diam5(i)=(1/36)*(-4*(-((f(i)^2*(-4+3*f(i))^3/(-5+4*f(i))^4)))^(1/2)+(f(i)*(68+9*f(i)*(-11+4*f(i))))/(5-4*f(i))^2);
    diam6(i)=-(20/9)*(-(f(i)^2*(-4+3*f(i))^3)/(-32+25*f(i))^4)^(1/2)+((f(i)*(-4+3*f(i))*(-104+75*f(i)))/(96-75*f(i))^2);
    diam2a(i)=(4/9)*((-8*(-(f(i)^2*(-4+3*f(i))^3)/(4+f(i))^4)^(1/2))+((f(i)*(80+9*(-8+f(i))*f(i)))/(4+f(i))^2));
end
figure
hold on
hold on
xlim([0.05,0.95]); % Auto-scale the x-axis
ylim([0.01,0.1]); % Auto-scale the y-axis
contourf(f,V,W')
colormap(flipud(bone))
colorbar

ylabel('V')
xlabel('f')

plot(f,diam2,'r','LineWidth',3)
plot(f,diam3,'y','LineWidth',3)
plot(f,diam4,'LineWidth',3)
plot(f,diam5,'LineWidth',3)
plot(f,diam6,'LineWidth',3)
annotation('textbox', [0.6, 0.72, 0.17, 0.07], 'String', 'Diameter 3', ...
    'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontSize', 12);
annotation('textbox', [0.6, 0.38, 0.17, 0.07], 'String', 'Diameter 4', ...
    'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontSize', 12);
lgd=legend('Diameter in the first period','Theoretical diameter = 2','Theoretical diameter = 3','Theoretical diameter = 4','Theoretical diameter = 5','Theoretical diameter = 6', 'Location', 'northwest');
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 18;
hold off

%export: 16 to 9 inches, font 16