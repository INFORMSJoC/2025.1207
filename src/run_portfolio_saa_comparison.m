% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

% Portfolio CVaR experiment.
clear all;
format long
clc;
namebiaohao=01;
filename=['portfolio_saa_comparison_',num2str(namebiaohao),'.xls'];
tic;
d=2;
phi=0.95;
miu=[1,1,1]';
SIGMA=[4,-2,-4;-2,5,0;-4,0,6];
theta0=[1,1,1];
alpha=1;calpha=1;
xunhuancishu=100;
colNum=1;

for c0=[1,2,4,10,20]
    rowNum=1;
syms y z u w;
x = sym('x', [1 d]);
ww = sym('ww', [1 d]);
w1=cos(x(1))^2;
w2=cos(x(2))^2*sin(x(1))^2;
w3=sin(x(2))^2*sin(x(1))^2;
w=[w1,w2,w3];
z=w*miu+sqrt(w*SIGMA*w.')*normpdf((sqrt(2)*erfinv(2*phi-1)-w*miu)/sqrt(w*SIGMA*w.'))/(1-phi);
% z=w*miu+sqrt(w*SIGMA*w.')/(sqrt(2*pi)*exp((erfinv(2*phi-1))^2)*(1-phi));
for i=1:d
    ww(i)=diff(z,x(i));
end
LB=[(9/400)*pi/2,(9/400)*pi/2];Delta=1-(9/400)*pi;
UB=[(1-9/400)*pi/2,(1-9/400)*pi/2];
digitsOld = digits(100);
thetatrue1=vpasolve([diff(z,x(1)) == 0,diff(z,x(2)) == 0], x,[0.7,1.1;0.7,1.1]);
thetatrue=zeros(1,d);
thetatrue=[double(thetatrue1.x1),double(thetatrue1.x2)];%subs(diff(z,x(1)),x,{0,0})
cvtrue=double(subs(z,x,thetatrue));
vtrue=double(subs(w*miu+sqrt(w*SIGMA*w.')*sqrt(2)*erfinv(2*phi-1),x,thetatrue));
cvjihe=zeros(c0,c0);
mincv=1000000;
for alphak=1:calpha
cvjihe=zeros(c0,c0);
for i=1:(c0)
for j=1:(c0)
% alpha=1.3;

% alpha=1+alphak/calpha;
thetajihe=[(9/400)*pi/2+i*Delta/c0,(9/400)*pi/2+j*Delta/c0];
alpha=vtrue;
% thetajihe=thetatrue;
N=200000;%1000,2300,4400,7300,10000,23000,44000,73000,100000
    rowNum=1;Xn=0;
    w1=cos(thetajihe(1))^2;
    w2=cos(thetajihe(2))^2*sin(thetajihe(1))^2;
    w3=sin(thetajihe(2))^2*sin(thetajihe(1))^2;
    w=[w1,w2,w3];
    for n=1:N
        Xn=Xn+max([w*miu+w*mvnrnd(zeros(d+1,d+1),SIGMA)*w'-alpha,0]);
    end
    cvjihe(i,j)=alpha+Xn/N/(1-phi);
    
end
end
[kk, linear_index] = min(cvjihe(:));
if kk<mincv
    mincv=kk;
    minalpha=alpha;
    [row, col] = ind2sub(size(cvjihe), linear_index);
    mintheta=[(9/400)*pi/2+row*Delta/c0,(9/400)*pi/2+col*Delta/c0];
end
toc;
end
disp(mincv);
disp(minalpha);
disp(mintheta);
disp((mincv-cvtrue)/cvtrue);
disp(calpha*c0*c0*N);


    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mincv,'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mintheta','sheet1',rowName);
    rowNum=rowNum+2;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,double((mincv-cvtrue)/cvtrue),'sheet1',rowName);
    rowNum=rowNum+1;
    colNum=colNum+1;
end
  %  1.427671443362778
  % 
  %  1.320000000000000
  % 
  %  0.825259957852980   0.871725666117692
  % 
  % -0.037494400978303


