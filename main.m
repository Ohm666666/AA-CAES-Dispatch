% 外层迭代控制主程序(子问题)
clear;clc;tic
%% 数据及参数
Data
Iteration.num = 0;      % 迭代次数
Iteration.time = [];    % 迭代用时
Iteration.RoCoF_u = []; % 每次迭代的RoCoF(功率突增扰动)
Iteration.RoCoF_d = []; % 每次迭代的RoCoF(功率突减扰动)
Iteration.FD60_u = [];  % 每次迭代的FD60(功率突增扰动)
Iteration.FD60_d = [];  % 每次迭代的FD60(功率突减扰动)
Iteration.MFD_u = [];   % 每次迭代的MFD(功率突增扰动)
Iteration.MFD_d = [];   % 每次迭代的MFD(功率突减扰动)
Iteration.cost = [];    % 每次迭代的目标函数
flag = zeros(T,2);      % FD60越限标志，越限时为1，否则为0
% parpool(20)
%% 优化割约束边界初值
f_u_SP = df60s*ones(T,1);   % FD60(功率突增扰动)边界
f_d_SP = df60s*ones(T,1);   % FD60(功率突减扰动)边界
%% 迭代
while 1
    MP_Dispatch(f_u_SP,f_d_SP)
    load MP_result.mat
    %% 功率突增扰动越限校验
    RoCoF_temp = zeros(T,1);    % 用于存放全时段RoCoF的临时变量
    FD60_temp = zeros(T,1);     % 用于存放全时段FD60的临时变量
    MFD_temp = zeros(T,1);      % 用于存放全时段MFD的临时变量
    parfor t = 1:T
        [RoCoF_temp(t),FD60_temp(t),MFD_temp(t)] = frequency(u_G_MP(t,:),P_G_MP(t,:),P_B_MP(t,:),P_L(t),step_W*sum(P_Wp(t,:))+step_L*P_L(t));
        if -FD60_temp(t)>df60s
            f_u_SP(t) = -f_u_MP(t)-max(-FD60_temp(t)-df60s,mss);
            flag(t,1) = 1;
        end
    end
    Iteration.RoCoF_u = [Iteration.RoCoF_u RoCoF_temp*fn];
    Iteration.FD60_u = [Iteration.FD60_u FD60_temp*fn];
    Iteration.MFD_u = [Iteration.MFD_u MFD_temp*fn];
    %% 功率突减扰动越限校验
    RoCoF_temp = zeros(T,1);    % 用于存放全时段RoCoF的临时变量
    FD60_temp = zeros(T,1);     % 用于存放全时段FD60的临时变量
    MFD_temp = zeros(T,1);      % 用于存放全时段MFD的临时变量
    parfor t = 1:T
        [RoCoF_temp(t),FD60_temp(t),MFD_temp(t)] = frequency(u_G_MP(t,:),P_G_MP(t,:),P_B_MP(t,:),P_L(t),-step_W*sum(P_Wp(t,:))-step_L*P_L(t));
        if FD60_temp(t)>df60s
            f_d_SP(t) = f_d_MP(t)-max(FD60_temp(t)-df60s,mss);
            flag(t,2) = 1;
        end
    end
    Iteration.RoCoF_d = [Iteration.RoCoF_d RoCoF_temp*fn];
    Iteration.FD60_d = [Iteration.FD60_d FD60_temp*fn];
    Iteration.MFD_d = [Iteration.MFD_d MFD_temp*fn];
    %% 迭代信息更新
    Iteration.num = Iteration.num+1;
    Iteration.time = [Iteration.time toc];
    Iteration.cost = [Iteration.cost [C_G_run_MP;C_G_ud_MP;C_W_MP]];
    disp(['已完成第',num2str(Iteration.num),'次迭代，当前总用时',num2str(Iteration.time(end)),'秒。'])
    %% 迭代结束判断
    if sum(flag,'all') == 0
        break
    else
        flag = zeros(T,2);  % 越限标志复位
    end
end
%% 结果保存
writematrix(u_G_MP,'result.xlsx','Sheet','G_result','range','B3:I98')
writematrix(P_G_MP,'result.xlsx','Sheet','G_result','range','J3:Q98')
writematrix(P_W_MP,'result.xlsx','Sheet','W_result','range','D3:E98')
writematrix(P_B_MP,'result.xlsx','Sheet','B_result','range','B3:C98')
writematrix(E_B_MP,'result.xlsx','Sheet','B_result','range','D3:E98')
writematrix(pha_MP,'result.xlsx','Sheet','Sys_result','range','B3:AN98')
writematrix(max(abs(P_tl_MP),[],2),'result.xlsx','Sheet','Sys_result','range','AO3:AO98')
writematrix(Iteration.RoCoF_u(:,end),'result.xlsx','Sheet','f_result','range','B3:B98')
writematrix(Iteration.FD60_u(:,end),'result.xlsx','Sheet','f_result','range','C3:C98')
writematrix(Iteration.MFD_u(:,end),'result.xlsx','Sheet','f_result','range','D3:D98')
writematrix(Iteration.RoCoF_d(:,end),'result.xlsx','Sheet','f_result','range','E3:E98')
writematrix(Iteration.FD60_d(:,end),'result.xlsx','Sheet','f_result','range','F3:F98')
writematrix(Iteration.MFD_d(:,end),'result.xlsx','Sheet','f_result','range','G3:G98')
writematrix(C_G_run_MP,'result.xlsx','Sheet','C_result','range','A2:A2')
writematrix(C_G_ud_MP,'result.xlsx','Sheet','C_result','range','B2:B2')
writematrix(C_W_MP,'result.xlsx','Sheet','C_result','range','C2:C2')
writematrix((1:Iteration.num)','result.xlsx','Sheet','Iteration','range','A2')
writematrix(Iteration.time','result.xlsx','Sheet','Iteration','range','B2')
writematrix(Iteration.cost','result.xlsx','Sheet','Iteration','range','C2')
save Iteration_result.mat Iteration