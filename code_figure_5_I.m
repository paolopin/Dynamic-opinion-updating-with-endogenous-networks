load('results_fig5-6.mat') % Comment this line to load different data

gg=grid_fine;
R=linspace(0.1,0.9,grid_fine);
for i=1:gg
    %Y1(i)=(f(i)*(4-3*f(i))*(n-3)*(n-1))/(36*n^2);
    %Y2(i)=(f(i)*(4-3*f(i))*(n-2)*(3*n-2))/(48*n^2);
    overlap(i) = (4*f(i)-3*f(i)^2)/16;
    %overlap(i) = (1/9)*(-8 * (-((f(i)^2 *(-4 + 3* f(i))^3)/((4 + f(i))^4)))^(1/2) + ( f(i) * (80 + 9 *(-8 + f(i)) *f(i)))/((4 + f(i) )^2) );
    overlap2(i) = f(i)/16;
    overlap3(i) = f(i) * ( 8 + 3*n)/(48*n);
    diam2(i)=(32*f(i)-36*f(i)^2+9*f(i)^3)/(9*f(i)^2)-(4/9)*((-(-64*f(i)^2+144*f(i)^3-108*f(i)^4+27*f(i)^5)/(f(i)^4)))^(1/2);
    diam3(i)=(1/36)*((f(i)*(-5+3*f(i))*(-4+3*f(i)))/(-1+f(i))^2-2*(-(f(i)^2*(-4+3*f(i))^3)/(-1+f(i))^4)^(1/2));
    diam4(i)=(f(i)*(-4+3*f(i))*(-40+27*f(i)))/(32-27*f(i))^2-12*(-(f(i)^2*(-4+3*f(i))^3)/(-32+27*f(i))^4)^(1/2);
    diam5(i)=(1/36)*(-4*(-((f(i)^2*(-4+3*f(i))^3/(-5+4*f(i))^4)))^(1/2)+(f(i)*(68+9*f(i)*(-11+4*f(i))))/(5-4*f(i))^2);
    diam6(i)=-(20/9)*(-(f(i)^2*(-4+3*f(i))^3)/(-32+25*f(i))^4)^(1/2)+((f(i)*(-4+3*f(i))*(-104+75*f(i)))/(96-75*f(i))^2);
    diam2a(i)=(4/9)*((-8*(-(f(i)^2*(-4+3*f(i))^3)/(4+f(i))^4)^(1/2))+((f(i)*(80+9*(-8+f(i))*f(i)))/(4+f(i))^2));
end
figure
hold on
xlim([0.05,0.95]); % Auto-scale the x-axis
ylim([0.01,0.1]); % Auto-scale the y-axis
%contourf(f,V,W')
%colormap(flipud(bone))
%surf(f,V,grid,'FaceColor','r','FaceAlpha',.4,'EdgeAlpha',.4)
% surf(Y,'FaceColor','b','FaceAlpha',.4,'EdgeAlpha',.4)
%plot(f,Y2,'LineWidth',6)
%plot(f,Y1,'LineWidth',6)
plot(f,overlap,'--','Color',[0.25,0.25,0.25],'LineWidth',3,'DisplayName','Proposition 5')
plot(f,overlap2,'k','LineWidth',3,'DisplayName','Proposition 6')

% Define the upper limit for shading (e.g., max y-value or top of the plot)
upper_limit = max(ylim); % This gets the upper y-limit of the current plot
lower_limit = min(ylim); % This gets the lower y-limit of the current plot
% Fill the area above the curve
fill([f, fliplr(f)], [overlap, upper_limit*ones(size(f))], [0.9, 0.9, 0.9], 'FaceAlpha', 0.65, 'EdgeColor', 'none');
fill([f, fliplr(f)], [overlap, fliplr(overlap2)], [0.7, 0.7, 0.7], 'FaceAlpha', 0.65, 'EdgeColor', 'none');
fill([f, fliplr(f)], [overlap2, lower_limit*ones(size(f))], [0.5, 0.5, 0.5], 'FaceAlpha', 0.65, 'EdgeColor', 'none');

ylabel('V')
xlabel('f')
for i=1:gg
    for j=1:gg
%         if grid(i,j)<3 %grid(i,j)<2 || ( grid(i,j)<6 && V(j)>overlap(i)-.01)
%             plot(f(i),V(j),'ok')
%         end
        if Y(i,j)>100
            plot(f(i),V(j),'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.95])
        end
            
    end
end
f_1=0.5;
sim2a=0.015;
sim2b=0.02;
sim3a=0.0350916;
sim3b=0.0350917;
sim4a=0.069;
sim4b=0.07;
% plot(f,diam2,'r','LineWidth',2)
% plot(f,diam3,'y','LineWidth',3)

% scatter(f_1,sim2a,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim2b,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim3a,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim3b,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim4a,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim4b,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
%plot(f,overlap3,'LineWidth',3)
% plot(f,diam4,'LineWidth',2)
% plot(f,diam5,'LineWidth',2)
% plot(f,diam6,'LineWidth',2)
lgd=legend('Proposition 5','Proposition 6','Consensus','Residual','Persistent disagreement','Simulation with disagreement', 'Location', 'northwest');
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 18;
hold off

figure
hold on
xlim([0.05,0.95]); % Auto-scale the x-axis
ylim([0.01,0.1]); % Auto-scale the y-axis
contourf(f,V,W')
colormap(flipud(bone))
%surf(f,V,grid,'FaceColor','r','FaceAlpha',.4,'EdgeAlpha',.4)
% surf(Y,'FaceColor','b','FaceAlpha',.4,'EdgeAlpha',.4)
%plot(f,Y2,'LineWidth',6)
%plot(f,Y1,'LineWidth',6)
plot(f,overlap,'--','Color',[0.25,0.25,0.25],'LineWidth',3,'DisplayName','Proposition 5')
plot(f,overlap2,'k','LineWidth',3,'DisplayName','Proposition 6')

% Define the upper limit for shading (e.g., max y-value or top of the plot)
upper_limit = max(ylim); % This gets the upper y-limit of the current plot
lower_limit = min(ylim); % This gets the lower y-limit of the current plot
% Fill the area above the curve
fill([f, fliplr(f)], [overlap, upper_limit*ones(size(f))], [0.9, 0.9, 0.9], 'FaceAlpha', 0.65, 'EdgeColor', 'none');
fill([f, fliplr(f)], [overlap, fliplr(overlap2)], [0.7, 0.7, 0.7], 'FaceAlpha', 0.65, 'EdgeColor', 'none');
fill([f, fliplr(f)], [overlap2, lower_limit*ones(size(f))], [0.5, 0.5, 0.5], 'FaceAlpha', 0.65, 'EdgeColor', 'none');

ylabel('V')
xlabel('f')
for i=1:gg
    for j=1:gg
%         if grid(i,j)<3 %grid(i,j)<2 || ( grid(i,j)<6 && V(j)>overlap(i)-.01)
%             plot(f(i),V(j),'ok')
%         end
        if Y(i,j)>100
            plot(f(i),V(j),'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.95])
        end
            
    end
end
f_1=0.5;
sim2a=0.015;
sim2b=0.02;
sim3a=0.0350916;
sim3b=0.0350917;
sim4a=0.069;
sim4b=0.07;
% plot(f,diam2,'r','LineWidth',2)
% plot(f,diam3,'y','LineWidth',3)

% scatter(f_1,sim2a,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim2b,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim3a,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim3b,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim4a,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
% scatter(f_1,sim4b,150,'s','MarkerEdgeColor','b','MarkerFaceColor','r')
%plot(f,overlap3,'LineWidth',3)
% plot(f,diam4,'LineWidth',2)
% plot(f,diam5,'LineWidth',2)
% plot(f,diam6,'LineWidth',2)
lgd=legend('Diameter in the first period','Proposition 5','Proposition 6','Consensus','Residual','Persistent disagreement','Simulation with disagreement', 'Location', 'northwest');
lgd.Title.String = 'Legend';
lgd.Title.FontSize = 18;
hold off

