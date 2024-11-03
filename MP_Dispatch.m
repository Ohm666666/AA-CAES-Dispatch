% 内层日前调度子程序(主问题)
function MP_Dispatch(f_u_SP,f_d_SP)
%% 数据及参数导入
yalmip('clear')
load data.mat
%% 变量定义
% 火电变量
u_G = binvar(T,G.N,'full');     % 火电机组运行状态
u_hp = binvar(T,G.N,'full');    % 高压定压段状态
u_sp = binvar(T,G.N,'full');    % 滑压段状态
u_lp = binvar(T,G.N,'full');    % 低压定压段状态
P_G = sdpvar(T,G.N,'full');     % 火电机组有功出力
P_Gu = sdpvar(T,G.N,'full');    % 功率突增扰动下火电机组PFR功率变化量
P_Gd = sdpvar(T,G.N,'full');    % 功率突减扰动下火电机组PFR功率变化量
Mc_u = sdpvar(T,G.N,'full');    % 功率突增扰动下McCormick辅助变量
Mc_d = sdpvar(T,G.N,'full');    % 功率突减扰动下McCormick辅助变量
% 电池储能变量
u_B = binvar(T,B.N,'full');     % 电池储能状态，放电为1，充电为0
P_B = sdpvar(T,B.N,'full');     % 电池储能功率，放电为正，充电为负
E_B = sdpvar(T,B.N,'full');     % 电池储能剩余能量
P_Bu = sdpvar(T,B.N,'full');    % 功率突增扰动下电池储能PFR功率变化量
P_Bd = sdpvar(T,B.N,'full');    % 功率突减扰动下电池储能PFR功率变化量
% 风电变量
P_W = sdpvar(T,W.N,'full');     % 风电场有功出力
% 系统变量
pha = sdpvar(T,NL.N,'full');    % 节点相位
P_tl = sdpvar(T,N_tl,'full');   % 线路有功功率(标幺值，方向：小节点指向大节点)
f_u = sdpvar(T,1,'full');       % 功率突增扰动下系统稳态频差
f_d = sdpvar(T,1,'full');       % 功率突减扰动下系统稳态频差
%% 目标函数
C_G_run = 0.25*sum(sum(repmat(G.Cr1,T,1).*P_G+repmat(G.Cr2,T,1).*u_G));     % 火电机组运行成本
C_G_ud = sum(sum(repmat(G.Cud,T-1,1).*abs(u_G(2:T,:)-u_G(1:T-1,:))));       % 火电机组启停成本
C_W = 0.25*sum(sum(penalty*(P_Wp-P_W)));                                    % 弃风惩罚成本
obj = C_G_run+C_G_ud+C_W+100*sum(f_d)-100*sum(f_u);
%% 火电机组运行约束
constraints = [];
% 火电机组状态约束
constraints = [constraints,u_hp+u_sp+u_lp==u_G];
% 火电机组出力约束
constraints = [constraints,u_G.*repmat(G.Pmin,T,1)<=P_G<=u_G.*repmat(G.Pmax,T,1)];
constraints = [constraints,repmat(G.Pmax.*G.pr1+m,T,1)-(1-u_hp)*M<=P_G<=repmat(G.Pmax,T,1)+(1-u_hp)*M];
constraints = [constraints,repmat(G.Pmax.*G.pr2+m,T,1)-(1-u_sp)*M<=P_G<=repmat(G.Pmax.*G.pr1,T,1)+(1-u_sp)*M];
constraints = [constraints,repmat(G.Pmin,T,1)-(1-u_lp)*M<=P_G<=repmat(G.Pmax.*G.pr2,T,1)+(1-u_lp)*M];
% 火电机组爬坡率约束
constraints = [constraints,-dt*repmat(G.ramp,T-1,1)-(2-u_G(2:T,:)-u_G(1:T-1,:))*M<=...
    P_G(2:T,:)-P_G(1:T-1,:)<=dt*repmat(G.ramp,T-1,1)+(2-u_G(2:T,:)-u_G(1:T-1,:))*M];
% 火电机组启停时间约束
for n = 1:G.N
    for t = 1:T-4*G.Tud(n)
        constraints = [constraints,sum(u_G(t+1:t+4*G.Tud(n),n))>=4*G.Tud(n)*(u_G(t+1,n)-u_G(t,n))];
        constraints = [constraints,sum(1-u_G(t+1:t+4*G.Tud(n),n))>=4*G.Tud(n)*(u_G(t,n)-u_G(t+1,n))];
    end
    for t = T-4*G.Tud(n)+1:T-1
        constraints = [constraints,sum(u_G(t+1:T,n))>=(T-t)*(u_G(t+1,n)-u_G(t,n))];
        constraints = [constraints,sum(1-u_G(t+1:T,n))>=(T-t)*(u_G(t,n)-u_G(t+1,n))];
    end
end
%% 电池储能运行约束
% 电池储能功率约束
constraints = [constraints,-(1-u_B)*M<=P_B<=repmat(B.Pmax,T,1)+(1-u_B)*M];
constraints = [constraints,-repmat(B.Pmax,T,1)-u_B*M<=P_B<=u_B*M];
% 电池储能能量约束
constraints = [constraints,E_B(1,:)==B.Eini];
constraints = [constraints,repmat(B.Emin,T,1)<=E_B<=repmat(B.Emax,T,1)];
constraints = [constraints,0.25*P_B(1:T-1,:)/B.eff-(1-u_B(1:T-1,:))*M<=...
    E_B(1:T-1,:)-E_B(2:T,:)<=0.25*P_B(1:T-1,:)/B.eff+(1-u_B(1:T-1,:))*M];
constraints = [constraints,0.25*P_B(1:T-1,:)*B.eff-u_B(1:T-1,:)*M<=...
    E_B(1:T-1,:)-E_B(2:T,:)<=0.25*P_B(1:T-1,:)*B.eff+u_B(1:T-1,:)*M];
constraints = [constraints,B.Eini-(1-u_B(T,:))*M<=E_B(T,:)-0.25*P_B(T,:)/B.eff<=B.Emax+(1-u_B(T,:))*M];
constraints = [constraints,B.Eini-u_B(T,:)*M<=E_B(T,:)-0.25*P_B(T,:)*B.eff<=B.Emax+u_B(T,:)*M];
%% 风电场运行约束
constraints = [constraints,0<=P_W<=P_Wp];
%% 系统运行约束
% 直流潮流约束(有名值)
for n = 1:N_tl
    constraints = [constraints,P_tl(:,n)==Ur^2*(pha(:,hn_tl(n))-pha(:,en_tl(n)))/X_tl(n)];
end
% 平衡节点(1节点)约束
constraints = [constraints,pha(:,1)==0];
% 线路功率约束
constraints = [constraints,-P_tlmax<=P_tl<=P_tlmax];
% 节点功率平衡约束(流入为正，流出为负，有名值)
for n = 1:NL.N
    -P_L*NL.pr(n);              % 添加有功负荷
    if NL.detp(n)=="N"
    elseif NL.detp(n)=="G"
        ans+P_G(:,NL.den(n));   % 添加火电机组有功出力
    elseif NL.detp(n)=="W"
        ans+P_W(:,NL.den(n));   % 添加风电场有功出力
    elseif NL.detp(n)=="B"
        ans+P_B(:,NL.den(n));   % 添加电池储能有功出力
    end
    for i = 1:length(NL.tln{n})
        if NL.tltp{n}(i)==1
            ans-P_tl(:,NL.tln{n}(i));   % 添加线路有功功率(本节点为线路首节点)
        elseif NL.tltp{n}(i)==0
            ans+P_tl(:,NL.tln{n}(i));   % 添加线路有功功率(本节点为线路末节点)
        end
    end
    constraints = [constraints,ans==0];
end
% 系统备用约束
constraints = [constraints,sum(u_G.*repmat(G.Pmax,T,1)-P_G,2)+sum(repmat(B.Pmax,T,1)-P_B,2)>=sum(P_Wer,2)+P_Ler];
constraints = [constraints,sum(u_G.*repmat(G.Pmax,T,1)-P_G,2)+sum(4*B.eff*(E_B-repmat(B.Emin,T,1))-P_B,2)>=sum(P_Wer,2)+P_Ler];
constraints = [constraints,sum(P_G-u_G.*repmat(G.Pmin,T,1),2)+sum(P_B+repmat(B.Pmax,T,1),2)>=sum(P_Wer,2)+P_Ler];
constraints = [constraints,sum(P_G-u_G.*repmat(G.Pmin,T,1),2)+sum(P_B-4*(E_B-repmat(B.Emax,T,1))/B.eff,2)>=sum(P_Wer,2)+P_Ler];
%% RoCoF约束
constraints = [constraints,step_W*sum(P_Wp,2)+step_L*P_L<=2*vfmax*sum(u_G.*repmat(G.Pmax.*G.H,T,1),2)];
%% 功率突增扰动下的FD60约束
% 火电机组功率变化量约束
constraints = [constraints,P_Gu<=-repmat(G.Pmax./G.R*G.hpr,T,1).*repmat(f_u+G.dz,1,G.N)+(1-u_hp)*M];
constraints = [constraints,P_Gu<=-repmat(G.Pmax./G.R.*G.lpr,T,1).*repmat(f_u+G.dz,1,G.N)+(1-u_lp)*M];
constraints = [constraints,P_Gu<=Mc_u+(1-u_sp)*M];
X_u = repmat(G.kpr./G.R,T,1).*P_G;
X_u_max = repmat(G.Pmax./G.R*G.hpr,T,1);
X_u_min = repmat(G.Pmax./G.R.*G.lpr,T,1);
Y_u = repmat(-f_u-G.dz,1,G.N);
Y_u_max = df60s-G.dz;
Y_u_min = 0;
constraints = [constraints,Mc_u>=X_u_min.*Y_u+Y_u_min*X_u-X_u_min*Y_u_min-(1-u_sp)*M];
constraints = [constraints,Mc_u>=X_u_max.*Y_u+Y_u_max*X_u-X_u_max*Y_u_max-(1-u_sp)*M];
constraints = [constraints,Mc_u<=X_u_min.*Y_u+Y_u_max*X_u-X_u_min*Y_u_max+(1-u_sp)*M];
constraints = [constraints,Mc_u<=X_u_max.*Y_u+Y_u_min*X_u-X_u_max*Y_u_min+(1-u_sp)*M];
constraints = [constraints,-u_sp*M<=Mc_u<=u_sp*M];
constraints = [constraints,P_Gu<=repmat(G.al*G.Pmax,T,1)+(1-u_G)*M];
constraints = [constraints,P_Gu<=repmat(G.Pmax,T,1)-P_G+(1-u_G)*M];
constraints = [constraints,0<=P_Gu<=u_G*M];
% 电池储能功率变化量约束
constraints = [constraints,P_Bu<=-repmat(B.kd.*B.Pmax,T,1).*repmat(f_u+B.dz,1,B.N)];
constraints = [constraints,P_Bu<=repmat(B.al*B.Pmax,T,1)];
constraints = [constraints,P_Bu<=repmat(B.Pmax,T,1)-P_B];
constraints = [constraints,P_Bu>=0];
% FD60约束
constraints = [constraints,step_W*sum(P_Wp,2)+step_L*P_L+kL*P_L.*f_u-sum(P_Gu,2)-sum(P_Bu,2)==0];
constraints = [constraints,max(G.dz,B.dz)<=-f_u<=df60s];
%% 功率突减扰动下的FD60约束
% 火电机组功率变化量约束
constraints = [constraints,-P_Gd<=repmat(G.Pmax./G.R*G.hpr,T,1).*repmat(f_d-G.dz,1,G.N)+(1-u_hp)*M];
constraints = [constraints,-P_Gd<=repmat(G.Pmax./G.R.*G.lpr,T,1).*repmat(f_d-G.dz,1,G.N)+(1-u_lp)*M];
constraints = [constraints,-P_Gd<=Mc_d+(1-u_sp)*M];
X_d = repmat(G.kpr./G.R,T,1).*P_G;
X_d_max = repmat(G.Pmax./G.R*G.hpr,T,1);
X_d_min = repmat(G.Pmax./G.R.*G.lpr,T,1);
Y_d = repmat(f_d-G.dz,1,G.N);
Y_d_max = df60s-G.dz;
Y_d_min = 0;
constraints = [constraints,Mc_d>=X_d_min.*Y_d+Y_d_min*X_d-X_d_min*Y_d_min-(1-u_sp)*M];
constraints = [constraints,Mc_d>=X_d_max.*Y_d+Y_d_max*X_d-X_d_max*Y_d_max-(1-u_sp)*M];
constraints = [constraints,Mc_d<=X_d_min.*Y_d+Y_d_max*X_d-X_d_min*Y_d_max+(1-u_sp)*M];
constraints = [constraints,Mc_d<=X_d_max.*Y_d+Y_d_min*X_d-X_d_max*Y_d_min+(1-u_sp)*M];
constraints = [constraints,-u_sp*M<=Mc_d<=u_sp*M];
constraints = [constraints,-P_Gd<=repmat(G.al*G.Pmax,T,1)+(1-u_G)*M];
constraints = [constraints,-P_Gd<=P_G-repmat(G.Pmin,T,1)+(1-u_G)*M];
constraints = [constraints,-u_G*M<=P_Gd<=0];
% 电池储能功率变化量约束
constraints = [constraints,-P_Bd<=repmat(B.kd.*B.Pmax,T,1).*repmat(f_d-B.dz,1,B.N)];
constraints = [constraints,-P_Bd<=repmat(B.al*B.Pmax,T,1)];
constraints = [constraints,-P_Bd<=P_B+repmat(B.Pmax,T,1)];
constraints = [constraints,P_Bd<=0];
% FD60约束
constraints = [constraints,-step_W*sum(P_Wp,2)-step_L*P_L+kL*P_L.*f_d-sum(P_Gd,2)-sum(P_Bd,2)==0];
constraints = [constraints,max(G.dz,B.dz)<=f_d<=df60s];
%% SP校验优化割约束
constraints = [constraints,-f_u<=f_u_SP];
constraints = [constraints,f_d<=f_d_SP];
%% 求解
ops = sdpsettings('solver','gurobi','verbose',1);
% ops.gurobi.TimeLimit = 200;
% ops.gurobi.MIPGap = 0.001;
ops.gurobi.TuneTimeLimit = 0;
solveinfo = optimize(constraints,obj,ops)
%% 结果保存
u_G_MP = round(value(u_G));
P_G_MP = u_G_MP.*max(min(value(P_G),repmat(G.Pmax,T,1)),repmat(G.Pmin,T,1));
u_B_MP = round(value(u_B));
P_B_MP = max(min(value(P_B),repmat(B.Pmax,T,1)),repmat(-B.Pmax,T,1));
E_B_MP = value(E_B);
P_W_MP = value(P_W);
pha_MP = value(pha);
P_tl_MP = value(P_tl);
f_u_MP = value(f_u);
f_d_MP = value(f_d);
C_G_run_MP = value(C_G_run);
C_G_ud_MP = value(C_G_ud);
C_W_MP = value(C_W);
save MP_result.mat u_G_MP P_G_MP u_B_MP P_B_MP E_B_MP P_W_MP pha_MP P_tl_MP...
    f_u_MP f_d_MP C_G_run_MP C_G_ud_MP C_W_MP solveinfo
end