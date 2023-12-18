% min_x (0.5 ||Ax-b||_2^2  + lambda ||x||_1)  / (lambda ||x||_{top_k})

clc;clear all;close all;
addpath('util','data','solver');


rand('seed',1);
randn('seed',1);


id_data = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];
result = [];
for idat = 1:length(id_data)
    idat
    randn('seed',idat);
    G = getdata_ica(id_data(idat));
    y = generate_y_DCSparseOpt(G);
    x0 =  0.1*randn(size(G,2),1);
    data_m = size(G,1);
    lambda = 0.1 / data_m;
    
    k = 100;
    MaxTime = 60;
    accuracy = 1e-10;
    speak = 0;
    % min_x (0.5 ||Gx-y||_2^2  + lambda ||x||_1)  / (lambda ||x||_{top_k})
    HandleObj = @(x)(0.5* norm(G*x-y)^2 + lambda * norm(x,1)) / (lambda*topksum(x,k));

    [x1,his1,ts1] = SpaseRecovery_DPA2(G,y,x0,lambda,k,speak,accuracy,MaxTime);
    [x2,his2,ts2] = SpaseRecovery_PGSA2(G,y,x0,lambda,k,speak,accuracy,MaxTime);
    [x3,his3,ts3] = SpaseRecovery_QTA(G,y,x0,lambda,k,speak,accuracy,MaxTime);
    [x4,his4,ts4] = SpaseRecovery_ParametricCD(G,y,x0,lambda,k,speak,accuracy,MaxTime);
    
 
    One = [];
    One.his1 = his1;
    One.his2 = his2;
    One.his3 = his3;
    One.his4 = his4;
    One.ts1 = ts1;
    One.ts2 = ts2;
    One.ts3 = ts3;
    One.ts4 = ts4;
    result{idat} = One;
    save(mfilename,'result')
 
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


pause(1)

end







