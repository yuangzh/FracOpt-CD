function example_one_dim
clc;clear all;close all;


pcolor.red      = [255,66,14]/256;    pcolor.red2     = [255,113,31]/256;   pcolor.red3     = [236,42,45]/256;
pcolor.blue     = [0,149,201]/256;    pcolor.blue2    = [114,214,238]/256;  pcolor.blue3    = [0,183,240]/256;  pcolor.blue4    = [0,192,192]/256;   pcolor.blue5    = [65,146,228]/256;
pcolor.green    = [0,170,77]/256;     pcolor.green2   = [81,157,28]/256;    pcolor.green3   = [36,178,76]/256;  pcolor.green4   = [192,192,0]/256;   pcolor.green5   = [76,181,60]/256;
pcolor.purple   = [125,66,210] /256;   pcolor.purple2  =[143,76,178] /256;
pcolor.crimson  = [212,22,118]/256 ;   pcolor.crimson2 = [245,147,202]/256; pcolor.crimson3  = [192,52,148]/256;
pcolor.orange   = [254,181,89]/256;   pcolor.orange2  = [255,129,0]/256;    pcolor.orange3  = [219,130,1]/256;
pcolor.gray     = [128,128,127]/256;
pcolor.yellow   = [255,204,51]/256;   pcolor.yellow2  = [248,235,46]/256;
pcolor.pink     = [255,130,160]/256;
pcolor.black     = [30,30,30]/256;
pcolor.def      = [256,256,230]/256;


 
xs = [-3.5:0.001:2.1];
f = @(x)(x+2)^2 / (abs(3*x+2)+1);
for i=1:length(xs),
    fs(i) = f(xs(i));
end

plot(xs,fs,'LineWidth',5,'color',pcolor.blue3);hold on;

cri = [-2;-2/3;0];
plot(cri(1),f(cri(1)));hold on;
text(cri(1),f(cri(1)),'x=-2','FontSize',25); hold on;
plot(cri(1),f(cri(1)),'r.','MarkerSize',20); hold on;

plot(cri(2),f(cri(2)),'LineWidth',5);hold on;
text(cri(2),f(cri(2)),'x=-2/3','FontSize',25); hold on;
plot(cri(2),f(cri(2)),'r.','MarkerSize',20);hold on;

plot(cri(3),f(cri(3)),'LineWidth',5);hold on;
text(cri(3),f(cri(3)),'x=0','FontSize',25); hold on;
plot(cri(3),f(cri(3)),'r.','MarkerSize',20);hold on;


fprintf('\n');
set(gcf,'paperpositionmode','auto')
print(sprintf('%s.eps',mfilename),'-depsc2','-loose');
print(sprintf('%s.png',mfilename),'-dpng');

