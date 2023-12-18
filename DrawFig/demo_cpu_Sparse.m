clear all
clc;clear all;close all;

load('demo_cpu_SpaseRecovery_Topk');
for idata = 1:211
One = result{idata};

pcolor = loadcolor;
figure('color','w')
myplot = @loglog;
myplot(One.ts1,One.his1,'-ks','LineWidth',3,'MarkerSize',3,'color', pcolor.blue); hold on;
myplot(One.ts2,One.his2,'--ms','LineWidth',3,'MarkerSize',3,'color',  pcolor.green); hold on;
myplot(One.ts3,One.his3,'-ro','LineWidth',3,'MarkerSize',3,'color', pcolor.yellow); hold on;
myplot(One.ts4,One.his4,'-ro','LineWidth',3,'MarkerSize',3,'color', pcolor.red); hold on;

hleg=legend('DPA','PGSA','QTA','PCD') ;
set(hleg,'FontSize',15,'FontWeight','normal');
set(hleg,'Fontname','times new Roman');
set(hleg,'Location','NorthEast');
set(gca,'Fontsize', 17);
xlabel('Time (seconds)','FontSize',20)
ylabel('Objective','FontSize',20)

% set(gca,'xtick',[5 10 20 30 40 50 60 70 80 90]);
% set(gca,'XTickLabel',{'64','128','256','512','1024'})
% set(gcf,'position', [300 100 600 500]);

hiss = [One.his1(:);One.his2(:);One.his3(:);One.his4(:)];
tss = [One.ts1(:);One.ts2(:);One.ts3(:);One.ts4(:)];
axis([min(tss) max(tss) 0.6*min(hiss) max(hiss)])

fprintf('\n');
set(gcf,'paperpositionmode','auto')
print(sprintf('%s_%d.eps',mfilename,idata),'-depsc2','-loose');
print(sprintf('%s_%d.png',mfilename,idata),'-dpng');



end




