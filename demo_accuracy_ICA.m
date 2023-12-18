

% min_x x'x / ||Ax||_4^2
clc;clear all;close all;
addpath('util','data','solver');
rand('seed',1);
randn('seed',1);

id_data = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];


times = 8;
result = [];

for idat = 1:length(id_data)
    
    avg_his1 = [];
    avg_his2 = [];
    avg_his3 = [];
    
    for it = 1:times,
%         randn('seed',it);
        G = getdata_ica(id_data(idat));
        x0 =  randn(size(G,2),1);
        MaxTime = 60;
        accuracy = 1e-10;
        speak = 0;
 
        [x1,his1] = ICA_PGA(G,x0,speak,accuracy,MaxTime);
        [x2,his2] = ICA_PowerMethod(G,x0,speak,accuracy,MaxTime);
        [x3,his3] = ICA_FractionalCD(G,x0,speak,accuracy,MaxTime);
        
        avg_his1(it) = min(his1); avg_his2(it) = min(his2); avg_his3(it) = min(his3);
    end
    
    One = [];
    One.his1 = [mean(avg_his1);std(avg_his1)];
    One.his2 = [mean(avg_his2);std(avg_his2)];
    One.his3 = [mean(avg_his3);std(avg_his3)];
    result{idat} = One; save(mfilename,'result')
    
    fprintf('%s & ',GetDataStr2(idat));
    fprintf('%.3f $\\pm$ %.3f & ',mean(avg_his1),std(avg_his1));
    fprintf('%.3f $\\pm$ %.3f & ',mean(avg_his2),std(avg_his2));
    fprintf('%.3f $\\pm$ %.3f  ',mean(avg_his3),std(avg_his3));
    fprintf('\\\\');
    fprintf('\n');
    
end

 