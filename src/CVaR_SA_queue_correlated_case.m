% Copyright (c) 2026 Qixin Wang, Hao Cao, and Jian-Qiang Hu
% SPDX-License-Identifier: MIT

clear all;
clc;

namebiaohao=635;
filename=['CVaR_SA_queue_regen_',num2str(namebiaohao),'.xlsx'];


scriptTimer=tic;
d=4;
phi=0.95;
CIalpha=0.05;

v=[0.1,0.2,0.1,0.3]';
v=v(1:d);
cons=[0.1,0.01];
A=[10,2,1,2;2,9,2,4;1,2,8,0;2,4,0,20];
thetahead=[7,8,9,10]';
theta0=[1,1,1,1];
v0=1;
cv0=1;
cu0=1;                 % regenerative 版本中，cujihe 存 m_n，而不是原单样本二阶矩
Gv0=[1,1,1,1];
Gcv0=[1,1,1,1];

% 保留原代码中的 stepsize 指数
alpha=0.9999;
beta=0.998;
gamma=0.997;
gammau=0.997;
varphi=0.996;
eta=0.995;
varsigma=1/6;
lamta=0.1;               % arrival rate lambda
xunhuancishu=1;
colNum=1;

% Waiting-time benchmark
s_q=double(v'*thetahead);
Ainvv=A\v;
b_q=double(v'*Ainvv);
kappa=cons(1)/(2*cons(2));

rho_fun=@(s) lamta.*s./(1+lamta.*s);
dcvar_fun=@(s) 1+log(rho_fun(s)./(1-phi)) ...
    +1./(1+lamta.*s);

root_fun=@(s) s-s_q+kappa*b_q*dcvar_fun(s);

s_lower=(1-phi)/(lamta*phi)*(1+1e-8);
s_true=fzero(root_fun,[s_lower,s_q]);

thetatrue=double(thetahead ...
    -kappa*dcvar_fun(s_true)*Ainvv);

rho_true=rho_fun(s_true);
vartrue=s_true*log(rho_true/(1-phi));
cvtrue=s_true*(1+log(rho_true/(1-phi)));

Gcvtrue=double((1-phi)*dcvar_fun(s_true)*v');

% 保留原来的投影区域，不改变算法模拟路径
tail_const=1-log(1-phi);
theta_box_center=tail_const*kappa*Ainvv+thetahead;

UB=theta_box_center+5;
LB=theta_box_center-5;

maxxxx=10000;          % 单个 regenerative cycle 的最大顾客数保护上限


% 参数区：这里默认沿用原 queue 脚本/文中参数。需要实验调参时只改这一块。
N1=10000;
alphanc=15;
betanc=2000;
gammanc=3;
gammaunc=10;
varphinc=400;
etanc=30;
cnc=20;


use_feasible_barcn=3;

for N=[25000]%1000,2300,5000,10000,23000,
    nTimer=tic;
    disp(N);
    rowNum=1;

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
        % 使用 N+1 是因为第 1 行存初值，第 n 次更新写入 n+1 行。
        thetajihe=zeros(N+1,d);
        vjihe=zeros(N+1,1);
        cvjihe=zeros(N+1,1);
        cujihe=zeros(N+1,1);
        cvupperjihe=zeros(N+1,1);
        cvlowerjihe=zeros(N+1,1);
        Gvjihe=zeros(N+1,d);
        Gcvjihe=zeros(N+1,d);
        taujihe=zeros(N+1,1);
        sumtaujihe=zeros(N+1,1);
        sigmajihe=zeros(N+1,1);

        thetajihe(1,:)=theta0(1:d);
        vjihe(1)=v0;
        cvjihe(1)=cv0;
        cujihe(1)=cu0;
        cvupperjihe(1)=cv0;
        cvlowerjihe(1)=cv0;
        Gvjihe(1,:)=Gv0(1:d);
        Gcvjihe(1,:)=Gcv0(1:d);

        for n=1:N
            cn=(cnc/(n+N1)^varsigma);
            betan=(betanc/(n+N1)^beta);
            gamman=(gammanc/(n+N1)^gamma);
            gammaun=(gammaunc/(n+N1)^gammau);
            varphin=(varphinc/(n+N1)^varphi);
            etan=(etanc/(n+N1)^eta);
            alphan=(alphanc/(n+N1)^alpha);

            % 原代码写作 norm(Gvjihe(n))，这里按行向量修正为 norm(Gvjihe(n,:))。
            % 这不是算法口径改变，只是避免 MATLAB 线性索引导致只取单个元素。
            barcn=cn/max([1,norm(Gvjihe(n,:))/sqrt(d),norm(Gcvjihe(n,:))/sqrt(d)]);
            incn=(2*binornd(1,0.5,d,1)-1);
% 只保证 theta +/- barcn*incn 仍然满足 v'*theta > 0
wc_now = v' * thetajihe(n,:)';
if wc_now <= 0
    error('Current theta invalid: v''*theta = %.6g at n = %d.', wc_now, n);
end

denom = sum(abs(v .* incn));
barcn = min(barcn, 0.49 * wc_now / denom);
barcn = max(barcn, realmin);
            % 可选合法性保护：保证 plus/minus 两个扰动点都能定义合法 M/M/1 服务率。
            if use_feasible_barcn==1
                wc_now=v'*thetajihe(n,:)';
                denom=sum(abs(v.*incn));
                if denom>0
                    barcn_feasible=0.49*wc_now/denom;
                    if barcn_feasible>0
                        barcn=min(barcn,barcn_feasible);
                    end
                end
                barcn=max(barcn,1e-8);
            end

            % -----------------------------------------------------------------
            % 6.3.2 regenerative-cycle simulation block
            % 注意：Xncycle / Xnzcycle / Xnfcycle 中存放的是 waiting before service，
            %       与原代码 Xn=(maxeventclock-daolaishijian-sn) 保持同一含义。
            % -----------------------------------------------------------------
            [Xncycle,taun]=MM1_regen_cycle_waiting(thetajihe(n,:)',v,lamta,maxxxx);
            [Xnzcycle,taunz]=MM1_regen_cycle_waiting(thetajihe(n,:)'+barcn*incn,v,lamta,maxxxx);
            [Xnfcycle,taunf]=MM1_regen_cycle_waiting(thetajihe(n,:)'-barcn*incn,v,lamta,maxxxx);

            taujihe(n+1)=taun;
            sumtaujihe(n+1)=sumtaujihe(n)+taun;

            vnz=vjihe(n)+barcn*Gvjihe(n,:)*incn;
            vnf=vjihe(n)-barcn*Gvjihe(n,:)*incn;

            % VaR recursion
            vjihe(n+1)=vjihe(n)+betan*sum(phi-(Xncycle<=vjihe(n)));

            % CVaR recursion
            cycle_tail_reward=sum(Xncycle.*(Xncycle>=vjihe(n)))/(1-phi);
            cvjihe(n+1)=cvjihe(n)+gamman*(cycle_tail_reward-taun*cvjihe(n));

            % m_n recursion
            cycle_noise=cycle_tail_reward-taun*cvjihe(n);
            cujihe(n+1)=cujihe(n)+gammaun*(cycle_noise^2-cujihe(n));

            % modified CI
            mean_tau=sumtaujihe(n+1)/n;
            sigmajihe(n+1)=sqrt(abs(cujihe(n+1)))/max(mean_tau,eps);
            cvupperjihe(n+1)=cvjihe(n+1)+norminv(1-CIalpha/2,0,1)*sqrt(gamman)*sigmajihe(n+1);
            cvlowerjihe(n+1)=cvjihe(n+1)-norminv(1-CIalpha/2,0,1)*sqrt(gamman)*sigmajihe(n+1);

            % G_v recursion：按 6.3.2 double-sum 公式展开，不改成 cycle average。
            Ivz=sum(Xnzcycle<=vnz);
            Ivf=sum(Xnfcycle<=vnf);
            Gv_cycle=(taunf*Ivz-taunz*Ivf)*incn'/(2*barcn);
            Gvjihe(n+1,:)=Gvjihe(n,:)-varphin*Gv_cycle;

            % G_fc recursion：按 6.3.2 double-sum 公式展开，不改成 R/tau 平均。
            Rz=sum(Xnzcycle.*(Xnzcycle>=vnz));
            Rf=sum(Xnfcycle.*(Xnfcycle>=vnf));
            Gcv_cycle=(taunf*Rz-taunz*Rf)*incn'/(2*barcn)-taunz*taunf*Gcvjihe(n,:);
            Gcvjihe(n+1,:)=Gcvjihe(n,:)+etan*Gcv_cycle;

            % theta recursion：保留原目标函数梯度结构。
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

        finalIndex=N+1;
        thetaguji(k,:)=thetajihe(finalIndex,:);
        thetarate(k)=norm(thetajihe(finalIndex,:)-thetatrue(1:d)')^2;
cvguji(k)=cvjihe(finalIndex);

% cvkkk 只保留给 CIN，CIN 的计算不动
cvkkk=cvtrue;

% 当前 theta_n 代入等待时间 CVaR 公式，仅用于 pseudo RMSE
s_now=double(v'*thetajihe(finalIndex,:)');
rho_now=lamta*s_now/(1+lamta*s_now);

if rho_now > (1-phi)
    cv_theta_now=s_now*(1+log(rho_now/(1-phi)));
else
    cv_theta_now=rho_now*s_now/(1-phi);
end

cvrmse(k)=(cvjihe(finalIndex)-cv_theta_now)^2;
        Gcvguji(k,:)=Gcvjihe(finalIndex,:);
        Gcvrmse(k)=norm(Gcvjihe(finalIndex,:)-Gcvtrue)^2;

        sigmank=sigmajihe(finalIndex);
        CIguji(k,:)=[cvjihe(finalIndex)-norminv(1-CIalpha/2,0,1)*sqrt(gamman)*sigmank, ...
                     cvjihe(finalIndex)+norminv(1-CIalpha/2,0,1)*sqrt(gamman)*sigmank];
        CIN(k)=(cvjihe(finalIndex)-cvkkk)/max(sqrt(gamman)*sigmank,eps);
        runtime_seconds(k)=toc(repeatTimer);
        fprintf('[CVaR_SA_queue_correlated_case] N=%d repeat=%d/%d runtime_seconds=%.6f\n',N,k,xunhuancishu,runtime_seconds(k));
    end
    total_runtime_seconds=toc(nTimer);

    disp('thetatrue');
    disp(thetatrue');
    disp('cvtrue');
    disp(cvtrue);
    disp('Gcvtrue');
    disp(Gcvtrue);

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
    disp('rmseGcv');
    disp(sqrt(mean(Gcvrmse)));

    disp('meanCI');
    disp(mean(CIguji,1));
    disp('stdCI');
    disp(std(CIguji,1));
    disp(mean(abs(CIN)<=norminv(1-CIalpha/2,0,1)));

    % % % % % 保存结果到 xlsx
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,mean(thetaguji,1)','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(thetaguji,1,1)','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,thetarela','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,sqrt(mean(thetarate)),'sheet1',rowName);
    % rowNum=rowNum+1;
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
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,mean(Gcvguji,1)','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(Gcvguji,1,1)','sheet1',rowName);
    % rowNum=rowNum+d;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,sqrt(mean(Gcvrmse)),'sheet1',rowName);
    % rowNum=rowNum+1;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,mean(CIguji,1)','sheet1',rowName);
    % rowNum=rowNum+2;
    % rowName = [char(colNum-1+65),num2str(rowNum)];
    % xlswrite(filename,std(CIguji,1)','sheet1',rowName);
    % rowNum=rowNum+1;
    % rowNum=rowNum+1;
    % if colNum==1
    %     xlswrite(filename,{'N','Repeat_count','Total_seconds','Mean_seconds','Std_seconds'},'Runtime_seconds','A1');
    %     xlswrite(filename,{'N','Repeat','Runtime_seconds'},'Runtime_seconds','A10');
    %     xlswrite(filename,{'单位：秒；Total_seconds 为该 N 的100次循环及初始化耗时，不含 Excel 导出。'},'Runtime_seconds','A9');
    % end
    % xlswrite(filename,[N,xunhuancishu,total_runtime_seconds,mean(runtime_seconds),std(runtime_seconds,1)],'Runtime_seconds',['A',num2str(colNum+1)]);
    % xlswrite(filename,[repmat(N,xunhuancishu,1),(1:xunhuancishu)',runtime_seconds],'Runtime_seconds',['A',num2str(11+(colNum-1)*xunhuancishu)]);
    % colNum=colNum+1;


    % 画图
    thetarelatu=zeros(N+1,d);
    Gcvrelatu=zeros(1,N+1);
    for n=1:(N+1)
        thetarelatu(n,:)=(thetajihe(n,:)-thetatrue(1:d)')./thetatrue(1:d)';
        Gcvrelatu(n)=norm(Gcvjihe(n,:));
    end
    cvrelatu=(cvjihe-cvtrue)/cvtrue;
    cvupperrelatu=(cvupperjihe-cvtrue)/cvtrue;
    cvlowerrelatu=(cvlowerjihe-cvtrue)/cvtrue;

    figure;
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

    toc(scriptTimer);
end

% ========================================================================
% Local function: simulate one regenerative busy cycle of an M/M/1 queue
% Output Xcycle stores waiting before service only, not system/sojourn time.
% This matches the original queue code:
% Xn = maxeventclock - daolaishijian - sn.
% ========================================================================
function [Xcycle,taucycle]=MM1_regen_cycle_waiting(theta,v,lamta,maxxxx)
    theta=theta(:);
    v=v(:);
    wc_theta=v'*theta;

    if wc_theta<=0
        error('Invalid theta: v''*theta must be positive to define the service rate.');
    end

    miu=lamta+1/wc_theta;

    if miu<=lamta
        error('Unstable M/M/1 queue: service rate must be larger than arrival rate.');
    end

    Xcycle=zeros(maxxxx,1);

    % 一个 busy cycle 从系统为空后的第一个顾客到达到系统再次为空为止。
    % 第一个顾客看到空系统，waiting before service = 0。
    sn=-log(1-rand)/miu;
    workload=sn;
    taucycle=1;
    Xcycle(taucycle)=0;

    while taucycle<maxxxx
        an=-log(1-rand)/lamta;

        % 下一个顾客到达时，前面积压工作量先减少 inter-arrival time。
        workload=workload-an;

        % 若该顾客到达前系统已经清空，则当前 busy cycle 结束。
        % 这个到达属于下一个 regenerative cycle，不纳入当前 cycle。
        if workload<=0
            break;
        end

        % 当前顾客到达时，已有 workload 的剩余服务工作量，
        % 因此它的排队等待时间就是 workload，不包含自己的 service time。
        sn=-log(1-rand)/miu;
        taucycle=taucycle+1;
        Xcycle(taucycle)=workload;

        % 当前顾客加入队列后，系统剩余工作量增加自己的服务时间。
        workload=workload+sn;
    end

    Xcycle=Xcycle(1:taucycle);
end
