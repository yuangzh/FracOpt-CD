% min_x x'x / ||Ax||_4^2
clc;clear all;close all;
addpath('util','data','solver');
rand('seed',1);
randn('seed',1);

id_data = [11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44];

id_data = [12 14];

result = [];

for idat = 1:length(id_data)
    idat
    randn('seed',idat);
    G = getdata_ica(id_data(idat));

    x0 =  randn(size(G,2),1);
    MaxTime = 60;
    accuracy = 1e-10;
    speak = 0;
    
    [x1,his1,ts1] = ICA_PGSA(G,x0,speak,accuracy,MaxTime);
    [x2,his2,ts2] = ICA_PowerMethod(G,x0,speak,accuracy,MaxTime);
    [x3,his3,ts3] = ICA_FractionalCD(G,x0,speak,accuracy,MaxTime);

    
    One = [];
    One.his1 = his1;
    One.his2 = his2;
    One.his3 = his3;
    One.ts1 = ts1;
    One.ts2 = ts2;
    One.ts3 = ts3;
    result{idat} = One;
    save(mfilename,'result')
    
    
end
