% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

% Lognormal heavy-tail CVaR experiment.

clear all;
clc;
namebiaohao=3553;
filename=['lognormal_results_',num2str(namebiaohao),'.xls'];
tic;
d=2;
phi=0.95;
cv0=1;
cu0=1;
v0=1;
Gv0=[1,1,1];
Gcv0=[1,2,1];
CIalpha=0.05;
CIwidthScale=1/sqrt(2);
miu=1;
sig=1;
theta0=[miu,sig];
alpha=0;
beta=0.998;
gamma=0.997;
gammau=0.997;
varphi=0.996;
eta=0.995;
varsigma=1/6;
% beta=2/3+0.0001/3;
% gamma=2/3-0.0001;
% gammau=2/3-0.0001;
% varphi=2/3;
% eta=2/3-0.0002;
% varsigma=1/9;
xunhuancishu=1;
rowNum=1;
thetatrue=zeros(1,d);
UB=zeros(1,d);
LB=zeros(1,d);
for i=1:d
    thetatrue(i)=theta0(i);
    UB(i)=thetatrue(i);
    LB(i)=thetatrue(i);
end
cvtrue=1/2*exp(miu+sig^2/2)*(1+erf(sig/sqrt(2)-erfinv(2*phi-1)))/(1-phi);
syms x y z w1 w2;
z=1/2*exp(x+y^2/2)*(1+erf(y/sqrt(2)-erfinv(2*phi-1)));
w1=subs(diff(z,x),{x,y},{miu,sig});
w2=subs(diff(z,y),{x,y},{miu,sig});
Gcvtrue=[double(w1),double(w2)];
N1=10000;
colNum=1;
for N=[500000]%560,1000,4400,10000,44000,100000,500000
    disp(N);
    rowNum=1;
    vjihe=zeros(N,1);
    cvjihe=zeros(N,1);
    cujihe=zeros(N,1);
    cvupperjihe=zeros(N,1);
    cvlowerjihe=zeros(N,1);
    Gvjihe=zeros(N,d);
    Gcvjihe=zeros(N,d);
    vjihe(1)=v0;
    cvjihe(1)=cv0;
    cujihe(1)=cu0;
    Gvjihe(1,:)=Gv0(1:d);
    Gcvjihe(1,:)=Gcv0(1:d);
    cvguji=zeros(xunhuancishu,1);
    CIguji=zeros(xunhuancishu,2);
    CIN=zeros(xunhuancishu,1);
    cvrmse=zeros(xunhuancishu,1);
    Gcvguji=zeros(xunhuancishu,d);
    Gcvrmse=zeros(xunhuancishu,1);
        for k=1:xunhuancishu
            for n=1:N
            %     disp(n);
                cn=(10/(n+100*N1)^varsigma);
                betan=(1000/(n+N1)^beta);
                gamman=(8/n^gamma);
                gammaun=(10/n^gammau);
                varphin=(1000/(n+0.1*N1)^varphi);
                etan=(3/(n+10*N1)^eta);
%                 cn=(10/(n+100*N1)^varsigma);
%                 betan=(10000/(n+N1)^beta);
%                 gamman=(10/n^gamma);
%                 gammaun=(10/n^gammau);
%                 varphin=(15000/(n+N1)^varphi);
%                 etan=(2/(n+10*N1)^eta);
                barcn=cn/max([1,norm(Gvjihe(n))/sqrt(d),norm(Gcvjihe(n))/sqrt(d)]);
                incn=(2*binornd(1,0.5,d,1)-1);
                W=normrnd(thetatrue(1),thetatrue(2));
                Wz=normrnd(thetatrue(1)+barcn*incn(1),thetatrue(2)+barcn*incn(2));
                Wf=normrnd(thetatrue(1)-barcn*incn(1),thetatrue(2)-barcn*incn(2));
                Xn=exp(W);
                Xnz=exp(Wz);
                Xnf=exp(Wf);
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
            cvguji(k)=cvjihe(n);
            cvkkk=cvtrue;
            cvrmse(k)=(cvjihe(n)-cvkkk)^2;
            Gcvguji(k,:)=Gcvjihe(n,:);
            Gcvrmse(k,:)=norm(Gcvjihe(n,:)-Gcvtrue)^2;
            CIguji(k,:)=[cvjihe(n)-CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            CIN(k)=(cvjihe(n)-cvkkk)*(1-phi)/(CIwidthScale*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2)));
            end
        end
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
    Gcvrela=(mean(Gcvguji,1)-Gcvtrue)./Gcvtrue;
    disp('relaGcv');
    disp(Gcvrela);
    disp('rmseGcv');
    disp(sqrt(mean(Gcvrmse)));
    
    disp('meanCI');
    disp(mean(CIguji,1));
    disp('stdCI');
    disp(std(CIguji,1));

    % zz=0;
    % for i =1:length(CIN)
    %     if -norminv(1-CIalpha/2)<=real(CIN(i)) && norminv(1-CIalpha/2)>=real(CIN(i))
    %         zz=zz+1;
    %     end
    % end
    % disp('CI converge probability');
    % disp(zz/xunhuancishu);  
%     edges = [-4:0.2:4];
%     h = histogram(CIN,edges,'Normalization','pdf');
%     hold on;
%     x=-4:0.01:4;
%     y=normpdf(x);
%     plot(x,y,'-','color','#FF4500');
%     legend('our algorithm','normal');
%     rowName = ['A',num2str(1)];
%     xlswrite(['lognormal_ci_statistics_',num2str(namebiaohao),'.xls'],CIN,'sheet1',rowName);

    Gcvrelatu=zeros(N+1,2);
    for n=1:(N+1)
        Gcvrelatu(n,:)=(Gcvjihe(n,:)-Gcvtrue)./Gcvtrue;
    end
    cvrelatu=(cvjihe-cvtrue)/cvtrue;
    cvupperrelatu=(cvupperjihe-cvtrue)/cvtrue;
    cvlowerrelatu=(cvlowerjihe-cvtrue)/cvtrue;
    plot(1:(N+1),Gcvrelatu(:,1),'-','color','#7E2F8E');
    hold on;
    plot(1:(N+1),Gcvrelatu(:,2),'-','color','m');
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
    legend('G_{cv,n}(1)','G_{cv,n}(2)','cv_n');
    axis([0 (N+1) -1 1]); 
    xlabel('number of simulation evaluations');
    ylabel('relative error');
    yticks([-1 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 1]);
    yticklabels({'-100%','-50%','-20%','-10%','0%','10%','20%','50%','100%'});
%     plot(1:(N+1),Gcvrelatu,'-','color','g');
% 
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,mean(cvguji,1),'sheet1',rowName);
%     rowNum=rowNum+1;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,std(cvguji,1),'sheet1',rowName);
%     rowNum=rowNum+1;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,double(cvrela),'sheet1',rowName);
%     rowNum=rowNum+1;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,sqrt(mean(cvrmse)),'sheet1',rowName);
% 
%     rowNum=rowNum+1;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,mean(CIguji,1)','sheet1',rowName);
%     rowNum=rowNum+2;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,std(CIguji,1)','sheet1',rowName);
%     rowNum=rowNum+2;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,mean(Gcvguji,1)','sheet1',rowName);
%     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,std(Gcvguji,1,1)','sheet1',rowName);
%     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,Gcvrela','sheet1',rowName);
%     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,sqrt(mean(Gcvrmse)),'sheet1',rowName);
%     rowNum=rowNum+1;
%     colNum=colNum+1;
    toc;
end
