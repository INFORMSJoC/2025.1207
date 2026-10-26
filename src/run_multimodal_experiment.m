% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

% Multimodal CVaR experiment.

clear all;
clc;
% parpool(32);
namebiaohao=400.2;
filename=['multimodal_results_',num2str(namebiaohao),'.xls'];
fileCIN=['multimodal_ci_statistics_',num2str(namebiaohao),'.xls'];
tic;
d=2;
phi=0.9;
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
sheding2=[1,2,2,3;0,0,0,0;1,2,2,4;1,2,2,4];
sheding3=[0,0,0,0;0,1,3,4];
% sheding3=[3,5;2,4;3,5;1,5];
theta0=[2,2];
v0=1;
cv0=1;
cu0=1;
Gv0=[1,1,1,1];
Gcv0=[1,1,1,1];
% alpha=1;
% beta=2/3+0.0001/3;
% gamma=2/3-0.0001;
% gammau=2/3-0.0001;
% varphi=2/3;
% eta=2/3-0.0002;
% varsigma=1/9;
alpha=0.9999;
beta=3/4+0.0001/3;
gamma=3/4-0.0001;
gammau=3/4-0.0001;
varphi=2/3;
eta=3/4-0.0002;
varsigma=1/8;

alphanc=0.5;
betanc=50;
gammanc=0.2;
gammaunc=0.2;
varphinc=100;
etanc=5;
cnc=1;

xunhuancishu=100;
rowNum=1;
x = sym('x', [1 d]);
w = sym('w', [1 d]);
syms y z u;
u1=0;u2=0;u3=0;
for i=1:1
    u2=u2+(x(i)-sheding2(i,1))*(x(i)-sheding2(i,2))*(x(i)-sheding2(i,3))*(x(i)-sheding2(i,4))+3;
end

for i=2:2
    u3=u3+(x(i)-sheding3(i,1))*(x(i)-sheding3(i,2))*(x(i)-sheding3(i,3))*(x(i)-sheding3(i,4))+3;
end

if fenbumoshi==1
	z=(int((y)*exp(-y),-log(1-phi),inf)*u2+int(exp(-y),-log(1-phi),inf)*u3)/(1-phi);
elseif fenbumoshi==2
    z=(int((y)*exp(-y^2/2),norminv(phi,0,1),inf)*u2+int(exp(-y^2/2),norminv(phi,0,1),inf)*u3)/(1-phi)/sqrt(2*pi);
end
for i=1:d
    w(i)=diff(z,x(i));
end
thetatrue=zeros(2,d);
UB=zeros(1,d);
LB=zeros(1,d);
for i=1:d
    thetatrue(1,i)=double(min(solve(diff(z,x(i)),x(i))));
    thetatrue(2,i)=double(max(solve(diff(z,x(i)),x(i))));
    UB(i)=max(thetatrue(:,i))+1;
    LB(i)=min(thetatrue(:,i))-1;
end
cvtrue1=double(subs(z,x,[thetatrue(1,1),thetatrue(1,2)]));
cvtrue2=double(subs(z,x,[thetatrue(1,1),thetatrue(2,2)]));
cvtrue3=double(subs(z,x,[thetatrue(2,1),thetatrue(1,2)]));
cvtrue4=double(subs(z,x,[thetatrue(2,1),thetatrue(2,2)]));
N1=1000000;
colNum=1;
for N=[450000,1000000]%1000,4500,10000,45000,100000,450000,1000000，100000,230000,450000,730000,1000000
    rowNum=1;
    thetajihe=zeros(N,d);
    vjihe=zeros(N,1);
    cvjihe=zeros(N,1);
    CIjihe=zeros(N,2);
    cujihe=zeros(N,1);
    cvupperjihe=zeros(N,1);
    cvlowerjihe=zeros(N,1);
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
                gamman=(gammanc/(n+0.1*N1)^gamma);
                gammaun=(gammaunc/(n+N1)^gammau);
                varphin=(varphinc/(n+N1)^varphi);
                etan=(etanc/(n+N1)^eta);
                alphan=(alphanc/(n)^alpha);
                barcn=cn/max([1,norm(Gvjihe(n))/sqrt(d),norm(Gcvjihe(n))/sqrt(d)]);
                incn=(2*binornd(1,0.5,1,1)-1);
                incn2=(2*binornd(1,0.5,1,1)-1);
                zhishu=(-log(1-rand));
                zhishuz=(-log(1-rand));
                zhishuf=(-log(1-rand));
                Xn=zhishu*((thetajihe(n,1)-sheding2(1,1)).*(thetajihe(n,1)-sheding2(1,2)).*(thetajihe(n,1)-sheding2(1,3)).*(thetajihe(n,1)-sheding2(1,4))+3);
                Xn=Xn+(thetajihe(n,2)-sheding3(2,1)).*(thetajihe(n,2)-sheding3(2,2)).*(thetajihe(n,2)-sheding3(2,3)).*(thetajihe(n,2)-sheding3(2,4))+3;
                Xnz=zhishuz*sum((thetajihe(n,1)+barcn*incn2'-sheding2(1,1)').*(thetajihe(n,1)+barcn*incn2'-sheding2(1,2)').*(thetajihe(n,1)+barcn*incn'-sheding2(1,3)').*(thetajihe(n,1)+barcn*incn'-sheding2(1,4)')+3);
                Xnz=Xnz+(thetajihe(n,2)+barcn*incn'-sheding3(2,1)').*(thetajihe(n,2)+barcn*incn'-sheding3(2,2)').*(thetajihe(n,2)+barcn*incn2'-sheding3(2,3)').*(thetajihe(n,2)+barcn*incn2'-sheding3(2,4)')+3;                
                Xnf=zhishuf*sum((thetajihe(n,1)-barcn*incn'-sheding2(1,1)').*(thetajihe(n,1)-barcn*incn'-sheding2(1,2)').*(thetajihe(n,1)-barcn*incn'-sheding2(1,3)').*(thetajihe(n,1)-barcn*incn'-sheding2(1,4)')+3);
                Xnf=Xnf+(thetajihe(n,2)-barcn*incn2'-sheding3(2,1)').*(thetajihe(n,2)-barcn*incn2'-sheding3(2,2)').*(thetajihe(n,2)-barcn*incn2'-sheding3(2,3)').*(thetajihe(n,2)-barcn*incn2'-sheding3(2,4)')+3;                
                vnz=vjihe(n)+barcn*Gvjihe(n,:)*[incn,incn2]';
                vnf=vjihe(n)-barcn*Gvjihe(n,:)*[incn,incn2]';
                if Xn<=vjihe(n)
                    vjihe(n+1)=vjihe(n)+betan*(phi-1);
                else
                    vjihe(n+1)=vjihe(n)+betan*(phi);
                end
                if Xn>=vjihe(n)
                    % disp(Xn);
                    cvjihe(n+1)=cvjihe(n)+gamman*(Xn/(1-phi)-cvjihe(n));
                    cujihe(n+1)=cujihe(n)+gammaun*(Xn^2-cujihe(n));
                else
                    cvjihe(n+1)=cvjihe(n)+gamman*(-cvjihe(n));
                    cujihe(n+1)=cujihe(n)+gammaun*(-cujihe(n));
                end
                cvupperjihe(n+1)=cvjihe(n+1)+norminv(1-CIalpha/2,0,1)*sqrt(gamman*abs(cujihe(n+1)-(cvjihe(n+1)*(1-phi))^2))/(1-phi);
                cvlowerjihe(n+1)=cvjihe(n+1)-norminv(1-CIalpha/2,0,1)*sqrt(gamman*abs(cujihe(n+1)-(cvjihe(n+1)*(1-phi))^2))/(1-phi);
                if Xnz<=vnz
                    if Xnf<=vnf
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(1-1)*[incn,incn2]/(2*barcn);
                    else
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(1-0)*[incn,incn2]/(2*barcn);
                    end
                else
                    if Xnf<=vnf
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(0-1)*[incn,incn2]/(2*barcn);
                    else
                        Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*(0-0)*[incn,incn2]/(2*barcn);
                    end
                end
                if Xnz>=vnz
                    if Xnf>=vnf
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*1-Xnf*1)*[incn,incn2]/(2*barcn)-Gcvjihe(n,:));
                    else
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*1-Xnf*0)*[incn,incn2]/(2*barcn)-Gcvjihe(n,:));
                    end
                else
                    if Xnf>=vnf
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*0-Xnf*1)*[incn,incn2]/(2*barcn)-Gcvjihe(n,:));
                    else
                        Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*((Xnz*0-Xnf*0)*[incn,incn2]/(2*barcn)-Gcvjihe(n,:));
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
            end
            thetaguji(k,:)=thetajihe(n,:);
            thetarate(k)=min([norm(thetajihe(n,:)-[thetatrue(1,1),thetatrue(1,2)]),norm(thetajihe(n,:)-[thetatrue(1,1),thetatrue(2,2)]),norm(thetajihe(n,:)-[thetatrue(2,1),thetatrue(1,2)]),norm(thetajihe(n,:)-[thetatrue(1,2),thetatrue(2,2)])])^2;
            cvguji(k)=cvjihe(n);
            cvkkk=subs(z,x,thetajihe(n,:));
            cvrmse(k)=(cvjihe(n)-cvkkk)^2;
            Gcvguji(k,:)=Gcvjihe(n,:);
            Gcvrmse(k,:)=norm(Gcvjihe(n,:)-(1-phi)*subs(w,x,thetajihe(n,:)))^2;
            CIguji(k,:)=[cvjihe(n)-norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            CIN(k)=(cvjihe(n)-cvtrue1)/sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))*(1-phi);
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
            CIN(k)=(cvjihe(n)-cvtrue1)/sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))*(1-phi);
        end
    end
    gk1=0;gk2=0;gk3=0;gk4=0;thetaguji1=[0,0];thetaguji2=[0,0];thetaguji3=[0,0];thetaguji4=[0,0];cvguji1=[0];cvguji2=[0];cvguji3=[0];cvguji4=[0];
    for iii=1:xunhuancishu
        if norm(thetaguji(iii,1)-thetatrue(1,1))<norm(thetaguji(iii,1)-thetatrue(2,1)) && norm(thetaguji(iii,2)-thetatrue(1,2))<norm(thetaguji(iii,2)-thetatrue(2,2))
            gk1=gk1+1;
            thetaguji1(gk1,:)=thetaguji(iii,:);
            cvguji1(gk1)=cvguji(iii);
        elseif norm(thetaguji(iii,1)-thetatrue(1,1))<norm(thetaguji(iii,1)-thetatrue(2,1)) && norm(thetaguji(iii,2)-thetatrue(1,2))>=norm(thetaguji(iii,2)-thetatrue(2,2))
            gk2=gk2+1;
            thetaguji2(gk2,:)=thetaguji(iii,:);
            cvguji2(gk2)=cvguji(iii);
        elseif norm(thetaguji(iii,1)-thetatrue(1,1))>=norm(thetaguji(iii,1)-thetatrue(2,1)) && norm(thetaguji(iii,2)-thetatrue(1,2))<norm(thetaguji(iii,2)-thetatrue(2,2))
            gk3=gk3+1;
            thetaguji3(gk3,:)=thetaguji(iii,:);
            cvguji3(gk3)=cvguji(iii);
        elseif norm(thetaguji(iii,1)-thetatrue(1,1))>=norm(thetaguji(iii,1)-thetatrue(2,1)) && norm(thetaguji(iii,2)-thetatrue(1,2))>=norm(thetaguji(iii,2)-thetatrue(2,2))
            gk4=gk4+1;
            thetaguji4(gk4,:)=thetaguji(iii,:);
            cvguji4(gk4)=cvguji(iii);
        end
    end
    disp('meantheta');
    disp(mean(thetaguji1,1));
    disp(mean(thetaguji2,1));
    disp(mean(thetaguji3,1));
    disp(mean(thetaguji4,1));
    disp('stdtheta');
    disp(std(thetaguji1,1,1));
    disp(std(thetaguji2,1,1));
    disp(std(thetaguji3,1,1));
    disp(std(thetaguji4,1,1));
    thetarela1=(mean(thetaguji1,1)-[thetatrue(1,1),thetatrue(1,2)])./[thetatrue(1,1),thetatrue(1,2)];
    thetarela2=(mean(thetaguji2,1)-[thetatrue(1,1),thetatrue(2,2)])./[thetatrue(1,1),thetatrue(2,2)];
    thetarela3=(mean(thetaguji3,1)-[thetatrue(2,1),thetatrue(1,2)])./[thetatrue(2,1),thetatrue(1,2)];
    thetarela4=(mean(thetaguji4,1)-[thetatrue(2,1),thetatrue(2,2)])./[thetatrue(2,1),thetatrue(2,2)];
    disp('relatheta');
    disp(thetarela1);
    disp(thetarela2);
    disp(thetarela3);
    disp(thetarela4);
    disp('ratetheta');
    disp(sqrt(mean(thetarate)));
    disp('meancv');
    disp(mean(cvguji1));
    disp(mean(cvguji2));
    disp(mean(cvguji3));
    disp(mean(cvguji4));
    disp('stdcv');
    disp(std(cvguji1,1));
    disp(std(cvguji2,1));
    disp(std(cvguji3,1));
    disp(std(cvguji4,1));
    cvrela1=(mean(cvguji1)-cvtrue1)/cvtrue1;
    cvrela2=(mean(cvguji2)-cvtrue2)/cvtrue2;
    cvrela3=(mean(cvguji3)-cvtrue3)/cvtrue3;
    cvrela4=(mean(cvguji4)-cvtrue4)/cvtrue4;
    disp('relacv');
    disp(double(cvrela1));
    disp(double(cvrela2));
    disp(double(cvrela3));
    disp(double(cvrela4));
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
    
    zz=0;
    for i =1:length(CIN)
        if -norminv(1-CIalpha/2)<=real(CIN(i)) && norminv(1-CIalpha/2)>=real(CIN(i))
            zz=zz+1;
        end
    end
    disp('CI converge probability');
    disp(zz/xunhuancishu);
    edges = [-4:0.2:4];
    h = histogram(CIN,edges,'Normalization','pdf');
    hold on;
    x=-4:0.01:4;
    y=normpdf(x);
    plot(x,y,'-','color','#FF4500');
    legend('our algorithm','normal');
    rowName = ['A',num2str(1)];
    xlswrite(fileCIN,CIN,'sheet1',rowName);

    
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(thetaguji1,1)','sheet1',rowName);
    rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(thetaguji2,1)','sheet1',rowName);
    rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(thetaguji3,1)','sheet1',rowName);
    rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(thetaguji4,1)','sheet1',rowName);
    rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(thetaguji,1,1)','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,thetarela1','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,thetarela2','sheet1',rowName);
    % rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,sqrt(mean(thetarate)),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,mean(cvguji,1),'sheet1',rowName);
    rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(cvguji,1),'sheet1',rowName);
    % rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,double(cvrela1),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,double(cvrela2),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,double(cvrela3),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,double(cvrela4),'sheet1',rowName);
    rowNum=rowNum+1;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,sqrt(mean(cvrmse)),'sheet1',rowName);
    rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,mean(Gcvguji,1)','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(Gcvguji,1,1)','sheet1',rowName);
    % rowNum=rowNum+d;
%     rowName = ['A',num2str(rowNum)];
%     xlswrite(filename,num2str(Gcvrela),'sheet1',rowName);
%     rowNum=rowNum+d;
    rowName = [char(colNum-1+65),num2str(rowNum)];
    xlswrite(filename,sqrt(mean(Gcvrmse)),'sheet1',rowName);
    rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,mean(CIguji,1)','sheet1',rowName);
    % rowNum=rowNum+2;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(CIguji,1)','sheet1',rowName);
    % rowNum=rowNum+1;
    % rowNum=rowNum+1;
    colNum=colNum+1;


    toc;
end

% end

