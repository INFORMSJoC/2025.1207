% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

% Portfolio CVaR experiment.
clear all;
clc;
namebiaohao=360;
filename=['portfolio_results_',num2str(namebiaohao),'.xls'];
fileCIN=['portfolio_ci_statistics_',num2str(namebiaohao),'.xls'];
tic;
d=2;
phi=0.95;
CIalpha=0.05;
CIwidthScale=1/sqrt(2);
miu=[1,1,1]';
SIGMA=[4,-2,-4;-2,5,0;-4,0,6];
theta0=[1,1,1];
v0=1;
cv0=1;
cu0=1;
Gv0=[1,1,1,1];
Gcv0=[1,1,1,1];
alpha=0.9999;
beta=0.998;
gamma=0.997;
gammau=0.997;
varphi=0.996;
eta=0.995;
varsigma=1/6;

xunhuancishu=100;
colNum=1;

syms y z u w;
x = sym('x', [1 d]);
ww = sym('ww', [1 d]);
w1=cos(x(1))^2;
w2=cos(x(2))^2*sin(x(1))^2;
w3=sin(x(2))^2*sin(x(1))^2;
w=[w1,w2,w3];
z = w*miu + sqrt(w*SIGMA*w.')*normpdf(norminv(phi))/(1-phi);
% z=w*miu+sqrt(w*SIGMA*w.')/(sqrt(2*pi)*exp((erfinv(2*phi-1))^2)*(1-phi));
for i=1:d
    ww(i)=diff(z,x(i));
end
thetatrue=zeros(1,d);
theta_bar = 9*pi/800;
LB = theta_bar*ones(1,d);
UB = (pi/2-theta_bar)*ones(1,d);
digitsOld = digits(100);
thetatrue=vpasolve([diff(z,x(1)) == 0,diff(z,x(2)) == 0], x,[0.7, 1.1; 0.7, 1.1]);
thetatrue=double([vpa(thetatrue.x1),vpa(thetatrue.x2)]);%subs(diff(z,x(1)),x,{0,0})
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
gammanc=10;
gammaunc=10;
varphinc=1500;
etanc=500;
cnc=0.4;
N1=10000;
for N=[1000,2300,4400,7300,10000,23000,44000,73000,100000]%1000,2300,4400,7300,10000,23000,44000,73000,100000
    rowNum=1;
    thetajihe=zeros(N,d);
    vjihe=zeros(N,1);
    cvjihe=zeros(N,1);
    cujihe=zeros(N,1);
    cvupperjihe=zeros(N,1);
    cvlowerjihe=zeros(N,1);
    Gvjihe=zeros(N,d);
    Gcvjihe=zeros(N,d);
    CIjihe=zeros(N,2);
    thetajihe(1,:)=theta0(1:d);
    vjihe(1)=v0;
    cvjihe(1)=cv0;
    cujihe(1)=cu0;
    Gvjihe(1,:)=Gv0(1:d);
    Gcvjihe(1,:)=Gcv0(1:d);
    thetaguji=zeros(xunhuancishu,d);
    thetarate=zeros(xunhuancishu,1);
    cvguji=zeros(xunhuancishu,1);
    CIguji=zeros(xunhuancishu,2);
    CIN=zeros(xunhuancishu,1);
    cvrmse=zeros(xunhuancishu,1);
    Gcvguji=zeros(xunhuancishu,d);
    Gcvrmse=zeros(xunhuancishu,1);
        for k=1:xunhuancishu
            for n=1:N
            %     disp(n)
                cn=(cnc/(n+N1)^varsigma);
                betan=(betanc/(n+N1)^beta);
                gamman=(gammanc/(n+N1)^gamma);
                gammaun=(gammaunc/(n+N1)^gammau);
                varphin=(varphinc/(n+N1)^varphi);
                etan=(etanc/(n+N1)^eta);
                alphan=(alphanc/(n)^alpha);
                barcn=cn/max([1,norm(Gvjihe(n))/sqrt(d),norm(Gcvjihe(n))/sqrt(d)]);
                incn=(2*binornd(1,0.5,d,1)-1);
                w1=cos(thetajihe(n,1))^2;
                w2=cos(thetajihe(n,2))^2*sin(thetajihe(n,1))^2;
                w3=sin(thetajihe(n,2))^2*sin(thetajihe(n,1))^2;
                w=[w1,w2,w3];
                Xn = w*miu + w*mvnrnd(zeros(1,d+1),SIGMA,1)';
                w1=cos(thetajihe(n,1)-barcn*incn(1))^2;
                w2=cos(thetajihe(n,2)-barcn*incn(2))^2*sin(thetajihe(n,1)-barcn*incn(1))^2;
                w3=sin(thetajihe(n,2)-barcn*incn(2))^2*sin(thetajihe(n,1)-barcn*incn(1))^2;
                w=[w1,w2,w3];
                Xnf = w*miu + w*mvnrnd(zeros(1,d+1),SIGMA,1)';
                w1=cos(thetajihe(n,1)+barcn*incn(1))^2;
                w2=cos(thetajihe(n,2)+barcn*incn(2))^2*sin(thetajihe(n,1)+barcn*incn(1))^2;
                w3=sin(thetajihe(n,2)+barcn*incn(2))^2*sin(thetajihe(n,1)+barcn*incn(1))^2;
                w=[w1,w2,w3];
                Xnz = w*miu + w*mvnrnd(zeros(1,d+1),SIGMA,1)';
                vnz=vjihe(n)+barcn*Gvjihe(n,:)*incn;
                vnf=vjihe(n)-barcn*Gvjihe(n,:)*incn;
                if Xn<=vjihe(n)
                    vjihe(n+1)=vjihe(n)+betan*(phi-1);
                else
                    vjihe(n+1)=vjihe(n)+betan*(phi);
                end
                if Xn>=vjihe(n)
                    cvjihe(n+1)=cvjihe(n)+gamman*(Xn/(1-phi)-cvjihe(n));
                    cujihe(n+1)=cujihe(n)+gammaun*(Xn^2-cujihe(n));
                else
                    cvjihe(n+1)=cvjihe(n)+gamman*(-cvjihe(n));
                    cujihe(n+1)=cujihe(n)+gammaun*(-cujihe(n));
                end
                cvupperjihe(n+1)=cvjihe(n+1)+CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*abs(cujihe(n+1)-(cvjihe(n+1)*(1-phi))^2))/(1-phi);
                cvlowerjihe(n+1)=cvjihe(n+1)-CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*abs(cujihe(n+1)-(cvjihe(n+1)*(1-phi))^2))/(1-phi);
                if Xnz<=vnz
                    if Xnf<=vnf
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(1-1)*incn'/(2*barcn);
                    else
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(1-0)*incn'/(2*barcn);
                    end
                else
                    if Xnf<=vnf
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(0-1)*incn'/(2*barcn);
                    else
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(0-0)*incn'/(2*barcn);
                    end
                end
                if Xnz>=vnz
                    if Xnf>=vnf
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*1-Xnf*1)*incn'/(2*barcn)-Gcvjihe(n,:));
                    else
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*1-Xnf*0)*incn'/(2*barcn)-Gcvjihe(n,:));
                    end
                else
                    if Xnf>=vnf
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*0-Xnf*1)*incn'/(2*barcn)-Gcvjihe(n,:));
                    else
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*0-Xnf*0)*incn'/(2*barcn)-Gcvjihe(n,:));
                    end
                end
                thetan2=thetajihe(n,:)-alphan*(Gcvjihe(n,:)/(1-phi));
                for i = 1:d
                    if (thetan2(i)>=LB(i)) && (thetan2(i)<=UB(i))
                        thetajihe(n+1,i)=thetan2(i);
                    elseif (thetan2(i)<LB(i))
                        thetajihe(n+1,i)=LB(i);
                    else
                        thetajihe(n+1,i)=UB(i);
                    end
                end
                CIjihe(n,:)=[cvjihe(n)-CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            end
            thetaguji(k,:)=thetajihe(n,:);
            thetarate(k)=norm(thetajihe(n,:)-thetatrue(1:d))^2;
            cvguji(k)=cvjihe(n);
            cvkkk=subs(z,x,thetajihe(n,:));
            cvrmse(k)=(cvjihe(n)-cvkkk)^2;
            Gcvguji(k,:)=Gcvjihe(n,:);
            Gcvrmse(k,:)=norm(Gcvjihe(n,:)-(1-phi)*subs(ww,x,thetajihe(n,:)))^2;
            CIguji(k,:)=[cvjihe(n)-CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            CIN(k)=(cvjihe(n)-cvtrue)*(1-phi)/(CIwidthScale*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2)));
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

    disp('meancv');
    disp(mean(cvguji));
    disp('stdcv');
    disp(std(cvguji,1));
    cvrela=(mean(cvguji)-cvtrue)/cvtrue;
    disp('relacv');
    disp(double(cvrela));
    disp('rmsecv');
    disp(sqrt(mean(cvrmse)));

    disp('meanGcv');
    disp(mean(Gcvguji,1));
    disp('stdGcv');
    disp(std(Gcvguji,1,1));
    % Gcvrela=(mean(Gcvguji,1)-Gcvtrue)./Gcvtrue;
    % disp('relaGcv');
    % disp(Gcvrela);
    disp('rmseGcv');
    disp(sqrt(mean(Gcvrmse)));
    
    disp('meanCI');
    disp(mean(CIguji,1));
    disp('stdCI');
    disp(std(CIguji,1));

%     zz=0;
%     for i =1:length(CIN)
%         if -norminv(1-CIalpha/2)<=real(CIN(i)) && norminv(1-CIalpha/2)>=real(CIN(i))
%             zz=zz+1;
%         end
%     end
%     disp('CI converge probability');
%     disp(zz/xunhuancishu);
%     edges = [-4:0.2:4];
%     h = histogram(CIN,edges,'Normalization','pdf');
%     hold on;
%     x=-4:0.01:4;
%     y=normpdf(x);
%     plot(x,y,'-','color','#FF4500');
%     legend('our algorithm','normal');
%     rowName = ['A',num2str(1)];
%     xlswrite(fileCIN,CIN,'sheet1',rowName);
%     
    thetarelatu=zeros(N+1,d);
    for n=1:(N+1)
        thetarelatu(n,:)=(thetajihe(n,:)-thetatrue(1:d))./thetatrue(1:d);
    end
    cvrelatu=(cvjihe-cvtrue)/cvtrue;
    cvupperrelatu=(cvupperjihe-cvtrue)/cvtrue;
    cvlowerrelatu=(cvlowerjihe-cvtrue)/cvtrue;
    plot(1:(N+1),thetarelatu(:,1),'-','color','#32CD32');
    hold on;
    plot(1:(N+1),thetarelatu(:,2),'-','color','g');
    hold on;
    plot(1:(N+1),cvrelatu,'-','color','#FF8C00');
    hold on;
    plot(1:(N+1),0.1*ones(1,N+1),'--','color','#8B8B00');
    hold on;
    plot(1:(N+1),0*ones(1,N+1),'--','color','#8B8B00');
    hold on;
    plot(1:(N+1),-0.1*ones(1,N+1),'--','color','#8B8B00');
    hold on;
    fill([1:(N+1), fliplr(1:(N+1))], [cvupperrelatu', fliplr(cvlowerrelatu')], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    legend('\theta_{n}(1)','\theta_{n}(2)','cv_n');
    axis([0 (N+1) -1 1]); 
    xlabel('number of simulation evaluations');
    ylabel('relative error');
    yticks([-1 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 1]);
    yticklabels({'-100%','-50%','-20%','-10%','0%','10%','20%','50%','100%'});
%     
    
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
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(cvguji,1),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,std(cvguji,1),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,double(cvrela),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,sqrt(mean(cvrmse)),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(CIguji,1)','sheet1',rowName);
    rowNum=rowNum+2;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,std(CIguji,1)','sheet1',rowName);
    rowNum=rowNum+1;
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(Gcvguji,1)','sheet1',rowName);
    rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,std(Gcvguji,1,1)','sheet1',rowName);
    rowNum=rowNum+d;
%     rowName = ['A',num2str(rowNum)];
%     xlswrite(filename,num2str(Gcvrela),'sheet1',rowName);
%     rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,sqrt(mean(Gcvrmse)),'sheet1',rowName);
    colNum=colNum+1;

    toc;
end

