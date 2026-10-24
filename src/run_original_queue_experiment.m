% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

% Queueing CVaR experiment.
clear all;
clc;
namebiaohao=354;
filename=['original_queue_results_',num2str(namebiaohao),'.xls'];
scriptTimer=tic;
d=4;
phi=0.95;
CIalpha=0.05;
CIwidthScale=1/sqrt(2);
fenbumoshi=1;%1:exp,2:normal
v=[0.1,0.2,0.1,0.3]';
v=v(1:d);
cons=[0.1,0.01];
% A=[1,2,3,4];
% A=A'*A;
% B=0.01*eye(4);
% A=A+B;
A=[10,2,1,2;2,9,2,4;1,2,8,0;2,4,0,20];
thetahead=[7,8,9,10]';
theta0=[1,1,1,1];
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

xunhuancishu=1;
colNum=1;
lamta=0.1;

s_q=double(v'*thetahead);
Ainvv=A\v;
b_q=double(v'*Ainvv);
kappa=cons(1)/(2*cons(2));

rho_fun=@(s) lamta.*s./(1+lamta.*s);
dcvar_fun=@(s) 1 ...
    + log(rho_fun(s)./(1-phi)) ...
    + 1./(1+lamta.*s);

root_fun=@(s) s-s_q+kappa*b_q*dcvar_fun(s);

s_lower=(1-phi)/(lamta*phi)*(1+1e-8);
s_true=fzero(root_fun,[s_lower,s_q]);

thetatrue=thetahead-kappa*dcvar_fun(s_true)*Ainvv;

rho_true=rho_fun(s_true);
vartrue=s_true*log(rho_true/(1-phi));
cvtrue=s_true*(1+log(rho_true/(1-phi)));
Gcvtrue=(1-phi)*dcvar_fun(s_true)*v';


UB=thetatrue+5;
LB=thetatrue-5;

maxxxx=1000;
alphanc=12;
betanc=6000;
gammanc=6;
gammaunc=6;
varphinc=13000;
etanc=200;
cnc=20;


N1=10000;
for N=[25000]%1000,2300,5000,10000,23000,25000
    nTimer=tic;
    disp(N);%100,180,240,500,730,1000,2400,4500,5000,7300,10000,2500,5000,7500,
    rowNum=1;
    thetajihe=zeros(N,d);
    vjihe=zeros(N,1);
    cvjihe=zeros(N,1);
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
    runtime_seconds=zeros(xunhuancishu,1);
        for k=1:xunhuancishu
            repeatTimer=tic;
            for n=1:N
            %     disp(n)
                cn=(cnc/(n+N1)^varsigma);
                betan=(betanc/(n+N1)^beta);
                gamman=(gammanc/(n+N1)^gamma);
                gammaun=(gammaunc/(n+N1)^gammau);
                varphin=(varphinc/(n+N1)^varphi);
                etan=(etanc/(n+N1)^eta);
                alphan=(alphanc/(n+N1)^alpha);
                barcn=cn/max([1,norm(Gvjihe(n))/sqrt(d),norm(Gcvjihe(n))/sqrt(d)]);
                incn=(2*binornd(1,0.5,d,1)-1);

%                 Xn=(-log(1-rand))*v'*thetajihe(n,:)';
%                 Xnz=(-log(1-rand))*v'*(thetajihe(n,:)'+barcn*incn);
%                 Xnf=(-log(1-rand))*v'*(thetajihe(n,:)'-barcn*incn);

                    lamta=0.1;
                    miu=lamta+1/(v'*thetajihe(n,:)');
                    vn=rand;
                    un2=rand;
                    an=-log(1-vn)/lamta;
                    sn=-log(1-un2)/miu;
                    daolaishijian=an;
                    maxeventclock=an+sn;
                    for i = 1:(maxxxx-1)
                        vn=rand;
                        un2=rand;
                        sn=-log(1-vn)/miu;
                        an=-log(1-un2)/lamta;
                        daolaishijian=an+daolaishijian;
                        if  daolaishijian>=maxeventclock
                            maxeventclock=daolaishijian+sn;
                            Xn=0;
                        else
                            maxeventclock=maxeventclock+sn;
                            Xn=(maxeventclock-daolaishijian-sn);
                        end
                    end
                    lamta=0.1;
                    miu=lamta+1/(v'*(thetajihe(n,:)'+barcn*incn));
                    vn=rand;
                    un2=rand;
                    an=-log(1-vn)/lamta;
                    sn=-log(1-un2)/miu;
                    daolaishijian=an;
                    maxeventclock=an+sn;
                    for i = 1:(maxxxx-1)
                        vn=rand;
                        un2=rand;
                        sn=-log(1-vn)/miu;
                        an=-log(1-un2)/lamta;
                        daolaishijian=an+daolaishijian;
                        if  daolaishijian>=maxeventclock
                            maxeventclock=daolaishijian+sn;
                            Xnz=0;
                        else
                            maxeventclock=maxeventclock+sn;
                            Xnz=(maxeventclock-daolaishijian-sn);
                        end
                    end
                    lamta=0.1;
                    miu=lamta+1/(v'*(thetajihe(n,:)'-barcn*incn));
                    vn=rand;
                    un2=rand;
                    an=-log(1-vn)/lamta;
                    sn=-log(1-un2)/miu;
                    daolaishijian=an;
                    maxeventclock=an+sn;
                    for i = 1:(maxxxx-1)
                        vn=rand;
                        un2=rand;
                        sn=-log(1-vn)/miu;
                        an=-log(1-un2)/lamta;
                        daolaishijian=an+daolaishijian;
                        if  daolaishijian>=maxeventclock
                            maxeventclock=daolaishijian+sn;
                            Xnf=0;
                        else
                            maxeventclock=maxeventclock+sn;
                            Xnf=(maxeventclock-daolaishijian-sn);
                        end
                    end

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
                thetan2=thetajihe(n,:)-alphan*(cons(1)*Gcvjihe(n,:)/(1-phi)+cons(2)*2*(thetajihe(n,:)-thetahead')*A);
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
            thetarate(k)=norm(thetajihe(n,:)-thetatrue(1:d)')^2;
            cvguji(k)=cvjihe(n);
            s_now=v'*thetajihe(n,:)';
            rho_now=lamta*s_now/(1+lamta*s_now);
            
            cvkkk=s_now*(1+log(rho_now/(1-phi)));
            cvrmse(k)=(cvjihe(n)-cvkkk)^2;
            Gcvguji(k,:)=Gcvjihe(n,:);
            s_now=v'*thetajihe(n,:)';
            rho_now=lamta*s_now/(1+lamta*s_now);
            dcvar_ds_now=1+log(rho_now/(1-phi))+1/(1+lamta*s_now);
            Gcv_target=(1-phi)*dcvar_ds_now*v';
            
            Gcvrmse(k)=norm(Gcvjihe(n,:)-Gcv_target)^2;
            CIguji(k,:)=[cvjihe(n)-CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi),cvjihe(n)+CIwidthScale*norminv(1-CIalpha/2,0,1)*sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2))/(1-phi)];
            CIN(k)=(cvjihe(n)-cvkkk)*(1-phi)/sqrt(gamman*(cujihe(n)-(cvjihe(n)*(1-phi))^2));
            runtime_seconds(k)=toc(repeatTimer);
            fprintf('[run_original_queue_experiment] N=%d repeat=%d/%d runtime_seconds=%.6f\n',N,k,xunhuancishu,runtime_seconds(k));
        end
    total_runtime_seconds=toc(nTimer);
    disp('meantheta');
    disp(mean(thetaguji,1));
    disp('stdtheta');
    disp(std(thetaguji,1,1));
    thetarela=(mean(thetaguji,1)-thetatrue(1:d)')./thetatrue(1:d)';
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

    disp('!!!!!!!!');
    disp(sum(sort(CIN) >= -CIwidthScale*norminv(1-CIalpha/2,0,1) & sort(CIN) <= CIwidthScale*norminv(1-CIalpha/2,0,1)))
  
%     edges = [-4:0.2:4];
%     h = histogram(CIN,edges,'Normalization','pdf');
%     hold on;
%     x=-4:0.01:4;
%     y=normpdf(x);
%     plot(x,y,'-','color','#FF4500');
%     legend('our algorithm','normal');
%     
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,mean(thetaguji,1)','sheet1',rowName);
%     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,std(thetaguji,1,1)','sheet1',rowName);
%     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,thetarela','sheet1',rowName);
%     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,sqrt(mean(thetarate)),'sheet1',rowName);
%     rowNum=rowNum+1;
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
%     rowNum=rowNum+1;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,mean(Gcvguji,1)','sheet1',rowName);
%     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,std(Gcvguji,1,1)','sheet1',rowName);
%     rowNum=rowNum+d;
% %     rowName = ['A',num2str(rowNum)];
% %     xlswrite(filename,num2str(Gcvrela),'sheet1',rowName);
% %     rowNum=rowNum+d;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,sqrt(mean(Gcvrmse)),'sheet1',rowName);
%     rowNum=rowNum+1;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,mean(CIguji,1)','sheet1',rowName);
%     rowNum=rowNum+2;
%     rowName = [char(colNum-1+65),num2str(rowNum)];
%     xlswrite(filename,std(CIguji,1)','sheet1',rowName);
%     rowNum=rowNum+1;
%     rowNum=rowNum+1;

%     if colNum==1
%         xlswrite(filename,{'N','Repeat_count','Total_seconds','Mean_seconds','Std_seconds'},'Runtime_seconds','A1');
%         xlswrite(filename,{'Units: seconds; Total_seconds includes all replications and initialization for this N, excluding Excel export.'},'Runtime_seconds','A9');
%         xlswrite(filename,{'N','Repeat','Runtime_seconds'},'Runtime_seconds','A10');
%     end
%     xlswrite(filename,[N,xunhuancishu,total_runtime_seconds,mean(runtime_seconds),std(runtime_seconds,1)], ...
%         'Runtime_seconds',['A',num2str(colNum+1)]);
%     xlswrite(filename,[repmat(N,xunhuancishu,1),(1:xunhuancishu)',runtime_seconds], ...
%         'Runtime_seconds',['A',num2str(11+(colNum-1)*xunhuancishu)]);
%     colNum=colNum+1;
% 
% 
    thetarelatu=zeros(N+1,d);
    Gcvrelatu=zeros(1,N+1);
    for n=1:(N+1)
        thetarelatu(n,:)=(thetajihe(n,:)-thetatrue(1:d)')./thetatrue(1:d)';
        Gcvrelatu(n)=norm(Gcvjihe(n,:));
    end
    cvrelatu=(cvjihe-cvtrue)/cvtrue;
    cvupperrelatu=(cvupperjihe-cvtrue)/cvtrue;
    cvlowerrelatu=(cvlowerjihe-cvtrue)/cvtrue;
    plot(1:(N+1),thetarelatu(:,1),'-','color','#32CD32');
    hold on;
    plot(1:(N+1),thetarelatu(:,2),'-','color','g');
    hold on;
    plot(1:(N+1),thetarelatu(:,3),'-','color','#228B22');
    hold on;
    plot(1:(N+1),thetarelatu(:,4),'-','color','#98FB98');
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
    legend('\theta_{n}(1)','\theta_{n}(2)','\theta_{n}(3)','\theta_{n}(4)','cv_n');
    axis([0 (N+1) -1 1]); 
    xlabel('number of outer iterations');
    ylabel('relative error');
    yticks([-1 -0.5 -0.2 -0.1 0 0.1 0.2 0.5 1]);
    yticklabels({'-100%','-50%','-20%','-10%','0%','10%','20%','50%','100%'});
    % % % % plot(1:(N+1),Gcvrelatu,'-','color','g');
    toc(scriptTimer);
end

