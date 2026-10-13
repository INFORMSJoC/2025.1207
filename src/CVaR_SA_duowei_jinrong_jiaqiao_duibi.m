% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

%jinrong
clear all;
format long
clc;
namebiaohao=0;
filename=['CVaR_SA_example_jinrong',num2str(namebiaohao),'_jiaqiao_duibi.xls'];
tic;
d=2;
phi=0.95;
miu=[1,1,1]';
SIGMA=[4,-2,-4;-2,5,0;-4,0,6];
theta0=[1,1,1];
v0=1;
cv0=1;
alpha=0.9999;
beta=4/7;
% gamma=0.997;
varsigma=1/7;
xunhuancishu=1;
colNum=1;

syms y z u w;
x = sym('x', [1 d]);
ww = sym('ww', [1 d]);
w1=cos(x(1))^2;
w2=cos(x(2))^2*sin(x(1))^2;
w3=sin(x(2))^2*sin(x(1))^2;
w=[w1,w2,w3];
z=w*miu+sqrt(w*SIGMA*w.')*normpdf((sqrt(2)*erfinv(2*phi-1)-w*miu)/sqrt(w*SIGMA*w.'))/(1-phi);
for i=1:d
    ww(i)=diff(z,x(i));
end
LB=[(9/400)*pi/2,(9/400)*pi/2];
UB=[(1-9/400)*pi/2,(1-9/400)*pi/2];
digitsOld = digits(100);
thetatrue1=vpasolve([diff(z,x(1)) == 0,diff(z,x(2)) == 0], x,[0.7,1.1]);
thetatrue=zeros(1,d);
thetatrue=[double(thetatrue1.x1),double(thetatrue1.x2)];%subs(diff(z,x(1)),x,{0,0})
% xs=0:0.01:pi/2;
% ys=0:0.01:pi/2;
% zs=zeros(length(xs),length(ys));
% for i=1:length(xs)
%     for j=1:length(ys)
%         zs(i,j)=subs(z,x,{xs(i),ys(j)});
%     end
% end
% surf(xs,ys,zs);
cvtrue=double(subs(z,x,thetatrue));
alphanc=1;
betanc=75;
% gammanc=10;
cnc=0.4;
N1=10000;
for N=[1000,2300,4400,7300,10000,23000,44000,73000,100000]%1000,2300,4400,7300,10000,23000,44000,73000,100000
    rowNum=1;
    thetajihe=zeros(N,d);
    vjihe=zeros(N,1);
    % cvjihe=zeros(N,1);
    thetajihe(1,:)=theta0(1:d);
    vjihe(1)=v0;
    % cvjihe(1)=cv0;
    thetaguji=zeros(xunhuancishu,d);
    thetarate=zeros(xunhuancishu,1);
    % cvguji=zeros(xunhuancishu,1);
    % cvrmse=zeros(xunhuancishu,1);
        for k=1:xunhuancishu
            for n=1:N
            %     disp(n)
                cn=(cnc/(n+N1)^varsigma);
                betan=(betanc/(n+N1)^beta);
                % gamman=(gammanc/(n+N1)^gamma);
                % gammaun=(gammaunc/(n+N1)^gammau);
                alphan=(alphanc/(n)^alpha);
                barcn=cn;
                incn=(2*binornd(1,0.5,d,1)-1);
                w1=cos(thetajihe(n,1))^2;
                w2=cos(thetajihe(n,2))^2*sin(thetajihe(n,1))^2;
                w3=sin(thetajihe(n,2))^2*sin(thetajihe(n,1))^2;
                w=[w1,w2,w3];
                Xn=w*miu+w*mvnrnd(zeros(d+1,d+1),SIGMA)*w';
                w1=cos(thetajihe(n,1)-barcn*incn(1))^2;
                w2=cos(thetajihe(n,2)-barcn*incn(2))^2*sin(thetajihe(n,1)-barcn*incn(1))^2;
                w3=sin(thetajihe(n,2)-barcn*incn(2))^2*sin(thetajihe(n,1)-barcn*incn(1))^2;
                w=[w1,w2,w3];
                Xnf=w*miu+w*mvnrnd(zeros(d+1,d+1),SIGMA)*w';
                w1=cos(thetajihe(n,1)+barcn*incn(1))^2;
                w2=cos(thetajihe(n,2)+barcn*incn(2))^2*sin(thetajihe(n,1)+barcn*incn(1))^2;
                w3=sin(thetajihe(n,2)+barcn*incn(2))^2*sin(thetajihe(n,1)+barcn*incn(1))^2;
                w=[w1,w2,w3];
                Xnz=w*miu+w*mvnrnd(zeros(d+1,d+1),SIGMA)*w';
                if Xn<=vjihe(n)
                    vjihe(n+1)=vjihe(n)+betan*(phi-1);
                else
                    vjihe(n+1)=vjihe(n)+betan*(phi);
                end
                % if Xn>=vjihe(n)
                %     cvjihe(n+1)=cvjihe(n)+gamman*(Xn/(1-phi)-cvjihe(n));
                % else
                %     cvjihe(n+1)=cvjihe(n)+gamman*(-cvjihe(n));
                % end
                thetan2=thetajihe(n,:)-alphan*(max(Xnz,vjihe(n))-max(Xnf,vjihe(n)))*incn'/(2*barcn)/(1-phi);
                for i = 1:d
                    if (thetan2(i)>=LB(i)) && (thetan2(i)<=UB(i))
                        thetajihe(n+1,i)=thetan2(i);
                    elseif (thetan2(i)<LB(i))
                        thetajihe(n+1,i)=LB(i);
                    else
                        thetajihe(n+1,i)=UB(i);
                    end
                end
            end
            thetaguji(k,:)=thetajihe(n,:);
            thetarate(k)=norm(thetajihe(n,:)-thetatrue(1:d))^2;
            % cvguji(k)=cvjihe(n);
            % cvkkk=subs(z,x,thetajihe(n,:));
            % cvrmse(k)=(cvjihe(n)-cvkkk)^2;
        end
    disp('meantheta');
    disp(mean(thetaguji,1));
    disp('stdtheta');
    disp(std(thetaguji,1,1));
    thetarela=(mean(thetaguji,1)-thetatrue(1:d))./thetatrue(1:d);
    disp('relatheta');
    disp(double(thetarela));
    disp('ratetheta');
    disp(sqrt(mean(thetarate)));

    % disp('meancv');
    % disp(mean(cvguji));
    % disp('stdcv');
    % disp(std(cvguji,1));
    % cvrela=(mean(cvguji)-cvtrue)/cvtrue;
    % disp('relacv');
    % disp(double(cvrela));
    % disp('rmsecv');
    % disp(sqrt(mean(cvrmse)));
%     thetarelatu=zeros(N+1,d);
%     for n=1:(N+1)
%         thetarelatu(n,:)=(thetajihe(n,:)-thetatrue(1:d))./thetatrue(1:d);
%     end
%     % cvrelatu=(cvjihe-cvtrue)/cvtrue;
% 
%     plot(1:(N+1),thetarelatu(:,1),'-','color','#32CD32');
%     hold on;
%     plot(1:(N+1),thetarelatu(:,2),'-','color','g');
%     hold on;
%     plot(1:(N+1),0.1*ones(1,N+1),'--','color','#8B8B00');
%     hold on;
%     plot(1:(N+1),0*ones(1,N+1),'--','color','#8B8B00');
%     hold on;
%     plot(1:(N+1),-0.1*ones(1,N+1),'--','color','#8B8B00');
%     hold on;
%     legend('\theta_{n}(1)','\theta_{n}(2)');
%     axis([0 (N+1) -1 1]); 
%     xlabel('number of simulation evaluations');
%     ylabel('relative error');
%     yticks([-1 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 1]);
%     yticklabels({'-100%','-50%','-20%','-10%','0%','10%','20%','50%','100%'});
% %     
    
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(thetaguji,1)','sheet1',rowName);
    rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,std(thetaguji,1,1)','sheet1',rowName);
    rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,thetarela','sheet1',rowName);
    rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,sqrt(mean(thetarate)),'sheet1',rowName);
    rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,mean(cvguji,1),'sheet1',rowName);
    % rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(cvguji,1),'sheet1',rowName);
    % rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,double(cvrela),'sheet1',rowName);
    % rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,sqrt(mean(cvrmse)),'sheet1',rowName);
    % rowNum=rowNum+1;
     colNum=colNum+1;
    toc;
end

