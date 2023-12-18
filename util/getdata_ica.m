function [A] = getdata_ica(iwhich)


switch iwhich
    
    case 11
        m = 1024; n = 1024;
        x = load('E2006_5000_10000');  x = x.x;
    case 12
        m = 1024; n = 2048;
        x = load('E2006_5000_10000');  x = x.x;
    case 13
        m = 1024; n = 1024;
        x = load('E2006_5000_10000');  x = x.x;
    case 14
        m = 2048; n = 1024;
        x = load('E2006_5000_10000');  x = x.x;
        
        
    case 21
        m = 1024; n = 1024;
x = load('news20'); x = double(x.x);
    case 22
        m = 1024; n = 2048;
x = load('news20'); x = double(x.x);
    case 23
        m = 1024; n = 1024;
x = load('news20'); x = double(x.x);
    case 24
        m = 2048; n = 1024;
x = load('news20'); x = double(x.x);
        
        
    case 31
        m = 1024; n = 1024;
x = load('sector_train'); x = double(x.x);
    case 32
        m = 1024; n = 2048;
x = load('sector_train'); x = double(x.x);
    case 33
        m = 1024; n = 1024;
x = load('sector_train'); x = double(x.x);
    case 34
        m = 2048; n = 1024;
x = load('sector_train'); x = double(x.x);
        
        
    case 41
        m = 1024; n = 1024;
x = load('TDT2'); x = double(x.fea);
    case 42
        m = 1024; n = 2048;
x = load('TDT2'); x = double(x.fea);
    case 43
        m = 1024; n = 1024;
x = load('TDT2'); x = double(x.fea);
    case 44
        m = 2048; n = 1024;
x = load('TDT2'); x = double(x.fea);
     
end

x = double(x)+1e-5*randn(size(x));
A = RandomSel(x,m,n);
A = full(A);
A = A./norm(A,'fro');



function X = RandomSel(x,m,n)
[data_m,data_n]=size(x);
sub_set_m = randperm(data_m,m);
sub_set_n = randperm(data_n,n);
X = x(sub_set_m,sub_set_n);



function [A] = scaleA(A)
[m,n] = size(A);
seq = randperm(m*n,round(0.2*m*n));
A(seq) = A(seq)*100;
