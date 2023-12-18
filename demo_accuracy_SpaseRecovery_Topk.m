% min_x (0.5 ||Ax-b||_2^2  + lambda ||x||_1)  / (lambda ||x||_{top_k})
clc;clear all;close all;
addpath('util','data','solver');


rand('seed',10);
randn('seed',10);





id_data = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];

times = 2;
result = [];
for idat = 1:length(id_data)
    
    avg_his1 = [];
    avg_his2 = [];
    avg_his3 = [];
    avg_his4 = [];
    
    for it = 1:times,
        randn('seed',it);
        
        G = getdata_ica(id_data(idat));
        y = generate_y_DCSparseOpt(G);
        x0 =  0.1*randn(size(G,2),1);
        data_m = size(G,1);
        lambda =  0.1 / data_m;
        
        k = 100;
        MaxTime = 60;
        accuracy = 1e-10;
        speak = 0;
        % min_x (0.5 ||Gx-y||_2^2  + lambda ||x||_1)  / (lambda ||x||_{top_k})
        HandleObj = @(x)(0.5* norm(G*x-y)^2 + lambda * norm(x,1)) / (lambda*topksum(x,k));
        [x1,his1] = SpaseRecovery_DPA2(G,y,x0,lambda,k,speak,accuracy,MaxTime);
        [x2,his2] = SpaseRecovery_PGSA2(G,y,x0,lambda,k,speak,accuracy,MaxTime);
        [x3,his3] = SpaseRecovery_QTA(G,y,x0,lambda,k,speak,accuracy,MaxTime);
        [x4,his4] = SpaseRecovery_ParametricCD(G,y,x0,lambda,k,speak,accuracy,MaxTime);
        avg_his1(it) = min(his1); avg_his2(it) = min(his2); avg_his3(it) = min(his3); avg_his4(it) = min(his4);
        
    end
    
    One = [];
    One.his1 = [mean(avg_his1);std(avg_his1)];
    One.his2 = [mean(avg_his2);std(avg_his2)];
    One.his3 = [mean(avg_his3);std(avg_his3)];
    One.his4 = [mean(avg_his4);std(avg_his4)];
    result{idat} = One; save(mfilename,'result')
    
    
    fprintf('%s & ',GetDataStr2(idat));
    fprintf('%.3f $\\pm$ %.3f & ',mean(avg_his1),std(avg_his1));
    fprintf('%.3f $\\pm$ %.3f & ',mean(avg_his2),std(avg_his2));
    fprintf('%.3f $\\pm$ %.3f & ',mean(avg_his3),std(avg_his3));
    fprintf('%.3f $\\pm$ %.3f  ', mean(avg_his4),std(avg_his4));
    fprintf('\\\\');
    fprintf('\n');
    
end








