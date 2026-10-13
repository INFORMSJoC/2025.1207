% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

%example5
%exp
clear all;
clc;
% parpool(32);
namebiaohao=348;
filename=['CVaR_SA_example2',num2str(namebiaohao),'-206.xls'];
fileCIN=['CVaR_SA_example2CIN',num2str(namebiaohao),'-206.xls'];
tic;
d=4;
phi=0.95;
CIalpha=0.05;
fenbumoshi=1;%1:exp,2:normal
sc1=[0,0,0,0];
sc2=[1,1,1,1];
sc3=[1,1,1,1];
% cons1=[2,2,2,6];
cons1=[0,0,0,0];
cons2=[2,0,0,0];
cons3=[1,0,0,0];
% cons3=[2,2,2,6];
sheding1=[0,0;0,0;0,0;0,0];
% sheding1=[1,3;2,4;3,5;1,5];
sheding2=[2,2;1,1;4,4;3,3];
sheding3=[1,1;2,2;3,3;4,4];
% sheding3=[3,5;2,4;3,5;1,5];
theta0=[1,2,3,4];
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

alphanc=0.5;
betanc=4000;
gammanc=20;
gammaunc=20;
varphinc=12000;
etanc=100;
cnc=2;
xunhuancishu=100;
zhifangtugeduan=20;
x = sym('x', [1 d]);
w = sym('w', [1 d]);
syms y z u;
u1=0;u2=0;u3=0;
for i=1:d
    u1=u1+sc1(i)*((x(i)-sheding1(i,1))*(x(i)-sheding1(i,2))+cons1(i));
end
for i=1:d
    u2=u2+sc2(i)*((x(i)-sheding2(i,1))*(x(i)-sheding2(i,2))+cons2(i));
end
for i=1:d
    u3=u3+sc3(i)*((x(i)-sheding3(i,1))*(x(i)-sheding3(i,2))+cons3(i));
end
if fenbumoshi==1
	z=(int((y)*exp(-y),-log(1-phi),inf)*u2+int(exp(-y),-log(1-phi),inf)*u3)/(1-phi);
elseif fenbumoshi==2
    z=(int((y)*exp(-y^2/2),norminv(phi,0,1),inf)*u2+int(exp(-y^2/2),norminv(phi,0,1),inf)*u3)/(1-phi)/sqrt(2*pi);
end
for i=1:d
    w(i)=diff(z,x(i));
end
thetatrue=zeros(1,d);
UB=zeros(1,d);
LB=zeros(1,d);
for i=1:d
    thetatrue(i)=double(solve(diff(z,x(i)),x(i)));
    UB(i)=thetatrue(i)+1;
    LB(i)=thetatrue(i)-1;
end
cvtrue=double(subs(z,x,thetatrue(1:d)));
% for gammanc=1
N1=1000000;
colNum=1;
for N=[1000,4500,10000,45000,100000,450000,1000000]%1000,4500,10000,45000,100000,450000,1000000
    rowNum=1;
    thetajihe=zeros(N,d);
    vjihe=zeros(N,1);
    cvjihe=zeros(N,1);
    CIjihe=zeros(N,2);
    cujihe=zeros(N,1);
    Gvjihe=zeros(N,d);
    Gcvjihe=zeros(N,d);
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
    if fenbumoshi==1
        for k=1:xunhuancishu
            for n=1:N
            %     disp(n);
                cn=(cnc/(n+N1)^varsigma);
                betan=(betanc/(n+N1)^beta);
                gamman=(gammanc/(n+N1)^gamma);
                gammaun=(gammaunc/(n+N1)^gammau);
                varphin=(varphinc/(n+N1)^varphi);
                etan=(etanc/(n+N1)^eta);
                alphan=(alphanc/n^alpha);
                barcn=cn/max([1,norm(Gvjihe(n))/sqrt(d),norm(Gcvjihe(n))/sqrt(d)]);
                incn=(2*binornd(1,0.5,d,1)-1);
                zhishu=(-log(1-rand));
                zhishuz=(-log(1-rand));
                zhishuf=(-log(1-rand));
                Xn=zhishu*sum(sc2(1:d).*((thetajihe(n,:)-sheding2(1:d,1)').*(thetajihe(n,:)-sheding2(1:d,2)')+cons2(1:d)));
                Xn=Xn+sum(sc3(1:d).*((thetajihe(n,:)-sheding3(1:d,1)').*(thetajihe(n,:)-sheding3(1:d,2)')+cons3(1:d)));
                Xnz=zhishuz*sum(sc2(1:d).*((thetajihe(n,:)+barcn*incn'-sheding2(1:d,1)').*(thetajihe(n,:)+barcn*incn'-sheding2(1:d,2)')+cons2(1:d)));
                Xnz=Xnz+sum(sc3(1:d).*((thetajihe(n,:)+barcn*incn'-sheding3(1:d,1)').*(thetajihe(n,:)+barcn*incn'-sheding3(1:d,2)')+cons3(1:d)));
                Xnf=zhishuf*sum(sc2(1:d).*((thetajihe(n,:)-barcn*incn'-sheding2(1:d,1)').*(thetajihe(n,:)-barcn*incn'-sheding2(1:d,2)')+cons2(1:d)));
                Xnf=Xnf+sum(sc3(1:d).*((thetajihe(n,:)-barcn*incn'-sheding3(1:d,1)').*(thetajihe(n,:)-barcn*incn'-sheding3(1:d,2)')+cons3(1:d)));
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
                thetan2=thetajihe(n,:)-alphan*Gcvjihe(n,:)/(1-phi);
                for i = 1:d
                    if (thetan2(i)>=LB(i)) && (thetan2(i)<=UB(i))
                        thetajihe(n+1,i)=thetan2(i);
                    elseif (thetan2(i)<LB(i))
                        thetajihe(n+1,i)=LB(i);
                    else
                        thetajihe(n+1,i)=UB(i);
                    end
                end
                CIjihe(n,:)=[cvjihe(n)-norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            end
            thetaguji(k,:)=thetajihe(n,:);
            thetarate(k)=norm(thetajihe(n,:)-thetatrue(1:d))^2;
            cvguji(k)=cvjihe(n);
            cvkkk=subs(z,x,thetajihe(n,:));
            cvrmse(k)=(cvjihe(n)-cvkkk)^2;
            Gcvguji(k,:)=Gcvjihe(n,:);
            Gcvrmse(k,:)=norm(Gcvjihe(n,:)-(1-phi)*subs(w,x,thetajihe(n,:)))^2;
            CIguji(k,:)=[cvjihe(n)-norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            CIN(k)=(cvjihe(n)-cvkkk)/sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))*(1-phi);
        end
    elseif fenbumoshi==2
        for k=1:xunhuancishu
            for n=1:N
            %     disp(n);
                cn=(cnc/(n+N1)^varsigma);
                betan=(betanc/(n+N1)^beta);
                gamman=(gammanc/(n+N1)^gamma);
                gammaun=(gammaunc/(n+N1)^gammau);
                varphin=(varphinc/(n+N1)^varphi);
                etan=(etanc/(n+N1)^eta);
                alphan=(alphanc/n^alpha);
                barcn=cn/max([1,norm(Gvjihe(n))/sqrt(d),norm(Gcvjihe(n))/sqrt(d)]);
                incn=(2*binornd(1,0.5,d,1)-1);                
                zhishu=normrnd(0,1);
                zhishuz=normrnd(0,1);
                zhishuf=normrnd(0,1);
                Xn=zhishu*sum(sc2(1:d).*((thetajihe(n,:)-sheding2(1:d,1)').*(thetajihe(n,:)-sheding2(1:d,2)')+cons2(1:d)));
                Xn=Xn+sum(sc3(1:d).*((thetajihe(n,:)-sheding3(1:d,1)').*(thetajihe(n,:)-sheding3(1:d,2)')+cons3(1:d)));
                Xnz=zhishuz*sum(sc2(1:d).*((thetajihe(n,:)+barcn*incn'-sheding2(1:d,1)').*(thetajihe(n,:)+barcn*incn'-sheding2(1:d,2)')+cons2(1:d)));
                Xnz=Xnz+sum(sc3(1:d).*((thetajihe(n,:)+barcn*incn'-sheding3(1:d,1)').*(thetajihe(n,:)+barcn*incn'-sheding3(1:d,2)')+cons3(1:d)));
                Xnf=zhishuf*sum(sc2(1:d).*((thetajihe(n,:)-barcn*incn'-sheding2(1:d,1)').*(thetajihe(n,:)-barcn*incn'-sheding2(1:d,2)')+cons2(1:d)));
                Xnf=Xnf+sum(sc3(1:d).*((thetajihe(n,:)-barcn*incn'-sheding3(1:d,1)').*(thetajihe(n,:)-barcn*incn'-sheding3(1:d,2)')+cons3(1:d)));
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
                thetan2=thetajihe(n,:)-alphan*Gcvjihe(n,:)/(1-phi);
                for i = 1:d
                    if (thetan2(i)>=LB(i)) && (thetan2(i)<=UB(i))
                        thetajihe(n+1,i)=thetan2(i);
                    elseif (thetan2(i)<LB(i))
                        thetajihe(n+1,i)=LB(i);
                    else
                        thetajihe(n+1,i)=UB(i);
                    end
                end
                CIjihe(n,:)=[cvjihe(n)-norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            end
            thetaguji(k,:)=thetajihe(n,:);
            thetarate(k)=norm(thetajihe(n,:)-thetatrue(1:d))^2;
            cvguji(k)=cvjihe(n);
            cvkkk=subs(z,x,thetajihe(n,:));
            cvrmse(k)=(cvjihe(n)-cvkkk)^2;
            Gcvguji(k,:)=Gcvjihe(n,:);
            Gcvrmse(k,:)=norm(Gcvjihe(n,:)-(1-phi)*subs(w,x,thetajihe(n,:)))^2;
            CIguji(k,:)=[cvjihe(n)-norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            CIN(k)=(cvjihe(n)-cvtrue)/sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))*(1-phi);
        end
    end
    disp('meantheta');
    disp(mean(thetaguji,1));
    disp('stdtheta');
    disp(std(thetaguji,1,1));
    thetarela=(mean(thetaguji,1)-thetatrue(1:d))./thetatrue(1:d);
    disp('relatheta');
    disp(thetarela);
    disp('ratetheta');
    disp(sqrt(mean(thetarate)));

    disp('meancv');
    disp(mean(cvguji));
    disp('stdcv');
    disp(std(cvguji,1));
    cvrela=(mean(cvguji)-cvtrue)/cvtrue;
    disp('relacv');
    disp('!!!!!!!!!!!!!!!!!!!!!');
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
%     thetarelatu=zeros(N+1,d);
%     for n=1:(N+1)
%         thetarelatu(n,:)=(thetajihe(n,:)-thetatrue(1:d))./thetatrue(1:d);
%     end
%     cvrelatu=(cvjihe-cvtrue)/cvtrue;
%     plot(1:(N+1),thetarelatu(:,1),'-','color','b');
%     hold on;
%     plot(1:(N+1),thetarelatu(:,2),'-','color','y');
%     hold on;
%     plot(1:(N+1),thetarelatu(:,3),'-','color','m');
%     hold on;
%     plot(1:(N+1),thetarelatu(:,4),'-','color','c');
%     hold on;
%     plot(1:(N+1),cvrelatu,'-','color','r');
%     hold on;
%     plot(1:(N+1),0.1*ones(1,N+1),'--','color','g');
%     hold on;
%     plot(1:(N+1),0*ones(1,N+1),'--','color','g');
%     hold on;
%     plot(1:(N+1),-0.1*ones(1,N+1),'--','color','g');
%     hold on;
%     legend('\theta_{n}(1)','\theta_{n}(2)','\theta_{n}(3)','\theta_{n}(4)','cv_n');
%     axis([0 (N+1) -1 1]); 
%     xlabel('number of simulation evaluations');
%     ylabel('relative error');
%     yticks([-1 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 1]);
%     yticklabels({'-100%','-50%','-20%','-10%','0%','10%','20%','50%','100%'});

    % nbins = zhifangtugeduan;
    % edges = [-2:0.2:2];
    % h = histogram(CIN,edges,'Normalization','pdf');
    % hold on;
    % x=-4:0.01:4;
    % y=normpdf(x);
    % plot(x,y,'-','color','#FF4500');

    
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
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(CIguji,1)','sheet1',rowName);
    rowNum=rowNum+2;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,std(CIguji,1)','sheet1',rowName);
    rowNum=rowNum+1;
    rowNum=rowNum+1;
    colNum=colNum+1;


    toc;
end

% end

