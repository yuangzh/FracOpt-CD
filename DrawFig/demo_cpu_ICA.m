clear all
clc;clear all;close all;

load('demo_cpu_ICA');
for idata = 1:2
One = result{idata};

pcolor = loadcolor;
figure('color','w')
myplot = @loglog;
myplot(One.ts1,One.his1,'-ks','LineWidth',3,'MarkerSize',3,'color', pcolor.blue); hold on;
myplot(One.ts2,One.his2,'--ms','LineWidth',3,'MarkerSize',3,'color',  pcolor.purple); hold on;
myplot(One.ts3,One.his3,'-ro','LineWidth',3,'MarkerSize',3,'color', pcolor.red); hold on;
% myplot(One.ts4,One.his4,'-.b+','LineWidth',4,'MarkerSize',3,'color', pcolor.yellow); hold on;
% myplot(One.ts5,One.his5,'-.k*','LineWidth',4,'MarkerSize',3,'color', pcolor.red ); hold on;

hleg=legend('DPA','Power Method','FCD') ;
set(hleg,'FontSize',15,'FontWeight','normal');
set(hleg,'Fontname','times new Roman');
set(hleg,'Location','NorthEast');
set(gca,'Fontsize', 17);
xlabel('Time (seconds)','FontSize',20)
ylabel('Objective','FontSize',20)

% set(gca,'xtick',[5 10 20 30 40 50 60 70 80 90]);
% set(gca,'XTickLabel',{'64','128','256','512','1024'})
% set(gcf,'position', [300 100 600 500]);

hiss = [One.his1(:);One.his2(:);One.his3(:)];
tss = [One.ts1(:);One.ts2(:);One.ts3(:)];
axis([min(tss) max(tss) 1*min(hiss) 0.3*max(hiss)])

fprintf('\n');
set(gcf,'paperpositionmode','auto')
print(sprintf('%s_%d.eps',mfilename,idata),'-depsc2','-loose');
print(sprintf('%s_%d.png',mfilename,idata),'-dpng');



end




