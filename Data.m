% 数据输入
%% 火电机组参数
G.N = 8;                                                                    % 机组数量
G.Pmax = readmatrix('parameters.xlsx','Sheet','unit','Range','F3:M3');      % 额定有功功率
G.Pmin = readmatrix('parameters.xlsx','Sheet','unit','Range','F4:M4');      % 有功功率下限
G.Cud = readmatrix('parameters.xlsx','Sheet','unit','Range','F5:M5');       % 单次启停成本
G.Tud = readmatrix('parameters.xlsx','Sheet','unit','Range','F6:M6');       % 最小启停机时间
G.ramp = readmatrix('parameters.xlsx','Sheet','unit','Range','F7:M7');      % 爬坡率
G.Cr1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F8:M8');       % 运行成本系数1
G.Cr2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F9:M9');       % 运行成本系数2
G.pr1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F10:M10');     % 高压定压段与滑压段的交界功率
G.hpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F11:F11');     % 高压定压值(所有机组相同)
G.pr2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F12:M12');     % 低压定压段与滑压段的交界功率
G.kpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F13:M13');     % 滑压段斜率
G.lpr = readmatrix('parameters.xlsx','Sheet','unit','Range','F14:M14');     % 低压定压值
G.pr3 = readmatrix('parameters.xlsx','Sheet','unit','Range','F15:F15');     % 再热定压段与滑压段的交界功率(所有机组相同)
G.Tch1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F16:M16');    % 高压定压段高压缸前汽室容积时间常数反比例系数
G.Tchs = readmatrix('parameters.xlsx','Sheet','unit','Range','F17:M17');    % 滑压段高压缸前汽室容积时间常数值
G.Tch2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F18:M18');    % 低压定压段高压缸前汽室容积时间常数反比例系数
G.Trhc = readmatrix('parameters.xlsx','Sheet','unit','Range','F19:M19');    % 再热定压段再热器容积时间常数反比例系数
G.Trhs = readmatrix('parameters.xlsx','Sheet','unit','Range','F20:M20');    % 再热滑压段再热器容积时间常数值
G.osn = readmatrix('parameters.xlsx','Sheet','unit','Range','F21:M21');     % 高压缸功率自然过调系数额定值
G.osk3 = readmatrix('parameters.xlsx','Sheet','unit','Range','F22:F22');    % 高压缸功率自然过调系数三次项系数(所有机组相同)
G.osk2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F23:F23');    % 高压缸功率自然过调系数二次项系数(所有机组相同)
G.osk1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F24:F24');    % 高压缸功率自然过调系数一次项系数(所有机组相同)
G.osk0 = readmatrix('parameters.xlsx','Sheet','unit','Range','F25:F25');    % 高压缸功率自然过调系数常数项(所有机组相同)
G.Fhpn = readmatrix('parameters.xlsx','Sheet','unit','Range','F26:M26');    % 高压缸输出功率占比额定值
G.Fhp4 = readmatrix('parameters.xlsx','Sheet','unit','Range','F27:F27');    % 高压缸输出功率占比四次项系数(所有机组相同)
G.Fhp3 = readmatrix('parameters.xlsx','Sheet','unit','Range','F28:F28');    % 高压缸输出功率占比三次项系数(所有机组相同)
G.Fhp2 = readmatrix('parameters.xlsx','Sheet','unit','Range','F29:F29');    % 高压缸输出功率占比二次项系数(所有机组相同)
G.Fhp1 = readmatrix('parameters.xlsx','Sheet','unit','Range','F30:F30');    % 高压缸输出功率占比一次项系数(所有机组相同)
G.Fhp0 = readmatrix('parameters.xlsx','Sheet','unit','Range','F31:F31');    % 高压缸输出功率占比常数项(所有机组相同)
G.R = readmatrix('parameters.xlsx','Sheet','unit','Range','F32:M32');       % 调差系数
G.Tse = readmatrix('parameters.xlsx','Sheet','unit','Range','F33:M33');     % 伺服机构时间常数
G.Td = readmatrix('parameters.xlsx','Sheet','unit','Range','F34:M34');      % 汽包容积时间常数
G.Tsh = readmatrix('parameters.xlsx','Sheet','unit','Range','F35:M35');     % 过热器容积时间常数
G.ksh = readmatrix('parameters.xlsx','Sheet','unit','Range','F36:M36');     % 过热器及主汽管道流量系数
G.al = readmatrix('parameters.xlsx','Sheet','unit','Range','F37:F37');      % 功率限幅比例(所有机组相同)
G.dz = readmatrix('parameters.xlsx','Sheet','unit','Range','F38:F38');      % 调频死区标幺值(所有机组相同)
G.H = readmatrix('parameters.xlsx','Sheet','unit','Range','F39:M39');       % 惯性时间常数
%% 电池储能
B.N = 2;                                                                    % 数量
B.Pmax = readmatrix('parameters.xlsx','Sheet','unit','Range','P6:Q6');      % 额定功率
B.Emax = readmatrix('parameters.xlsx','Sheet','unit','Range','P7:Q7');      % 额定容量
B.Emin = readmatrix('parameters.xlsx','Sheet','unit','Range','P8:Q8');      % 能量下限
B.Eini = readmatrix('parameters.xlsx','Sheet','unit','Range','P9:Q9');      % 初始能量
B.eff = readmatrix('parameters.xlsx','Sheet','unit','Range','P10:P10');     % 效率(所有机组相同)
B.kv = readmatrix('parameters.xlsx','Sheet','unit','Range','P11:Q11');      % 虚拟惯量系数
B.kd = readmatrix('parameters.xlsx','Sheet','unit','Range','P12:Q12');      % 频率下垂系数
B.Tc = readmatrix('parameters.xlsx','Sheet','unit','Range','P13:Q13');      % 变流器时间常数
B.dz = readmatrix('parameters.xlsx','Sheet','unit','Range','P14:P14');      % 调频死区标幺值(所有机组相同)
B.al = readmatrix('parameters.xlsx','Sheet','unit','Range','P15:P15');      % 功率限幅比例(所有机组相同)
%% 风电及负荷预测功率
W.N = 2;                                                                    % 风电场数量
P_Wp = readmatrix('parameters.xlsx','Sheet','wind','Range','B3:C98');       % 风电预测出力
P_Wer = readmatrix('parameters.xlsx','Sheet','wind','Range','D3:E98');      % 风电预测误差
P_L = readmatrix('parameters.xlsx','Sheet','load','Range','B2:B97');        % 负荷预测功率
P_Ler = readmatrix('parameters.xlsx','Sheet','load','Range','C2:C97');      % 负荷预测误差
%% 线路参数
N_tl = 46;                                                                  % 线路数量
hn_tl = readmatrix('parameters.xlsx','Sheet','line','Range','B2:B47');      % 线路首节点编号
en_tl = readmatrix('parameters.xlsx','Sheet','line','Range','C2:C47');      % 线路末节点编号
X_tl = readmatrix('parameters.xlsx','Sheet','line','Range','D2:D47');       % 线路电抗(有名值)
R_tl = readmatrix('parameters.xlsx','Sheet','line','Range','E2:E47');       % 线路电阻(有名值)
G_tl = real(1./(R_tl+X_tl.*1i));                                            % 线路电导(有名值)
B_tl = imag(1./(R_tl+X_tl.*1i));                                            % 线路电纳(有名值)
%% 节点参数
NL.N = 39;                                                                  % 节点数量
NL.pr = readmatrix('parameters.xlsx','Sheet','node','Range','B3:B41');      % 节点有功负荷比例
NL.detp = readmatrix('parameters.xlsx','Sheet','node','Range','C3:C41','OutputType','string');   % 节点所接设备类型
NL.den = readmatrix('parameters.xlsx','Sheet','node','Range','D3:D41');     % 节点所接设备编号
NL.tln = cell(NL.N,1);                                                      % 节点相连线路编号
NL.tltp = cell(NL.N,1);                                                     % 节点相连线路类型，节点为线路首节点时取1，为线路末节点时取0
for i = 1:NL.N
    for n = 1:N_tl
        if hn_tl(n)==i
            NL.tln{i} = [NL.tln{i},n];
            NL.tltp{i} = [NL.tltp{i},1];
        elseif en_tl(n)==i
            NL.tln{i} = [NL.tln{i},n];
            NL.tltp{i} = [NL.tltp{i},0];
        end
    end
end
%% 其他参数
dt = readmatrix('parameters.xlsx','Sheet','other','Range','B2:B2');         % 单位调度时长
T = readmatrix('parameters.xlsx','Sheet','other','Range','B3:B3');          % 日前调度时段数
kL = readmatrix('parameters.xlsx','Sheet','other','Range','B4:B4');         % 负荷频率响应系数
fn = readmatrix('parameters.xlsx','Sheet','other','Range','B5:B5');         % 额定频率
vfmax = readmatrix('parameters.xlsx','Sheet','other','Range','B6:B6');      % 初始频率变化率限值标幺值
df60s = readmatrix('parameters.xlsx','Sheet','other','Range','B7:B7');      % 60s频差限值标幺值
dfmax = readmatrix('parameters.xlsx','Sheet','other','Range','B8:B8');      % 最大频差限值标幺值
dt_ni = readmatrix('parameters.xlsx','Sheet','other','Range','B9:B9');      % 数值积分步长
T_ni = readmatrix('parameters.xlsx','Sheet','other','Range','B10:B10');     % 数值积分总时长
N_ni = readmatrix('parameters.xlsx','Sheet','other','Range','B11:B11');     % 数值积分时步数
P_tlmax = readmatrix('parameters.xlsx','Sheet','other','Range','B12:B12');  % 线路容量上限
Ur = readmatrix('parameters.xlsx','Sheet','other','Range','B13:B13');       % 电压等级
Sr = readmatrix('parameters.xlsx','Sheet','other','Range','B14:B14');       % 功率基准值
penalty = readmatrix('parameters.xlsx','Sheet','other','Range','B15:B15');  % 弃风惩罚成本系数
step_W = readmatrix('parameters.xlsx','Sheet','other','Range','B16:B16');   % 风电阶跃比例
step_L = readmatrix('parameters.xlsx','Sheet','other','Range','B17:B17');   % 负荷阶跃比例
mss = readmatrix('parameters.xlsx','Sheet','other','Range','B18:B18');      % 60s频差上限修正的最大步长(maximum step size)
M = 1000000;                                                                % 足够大的正数
m = 0.0001;                                                                 % 足够小的正数
%% 结果输出
clear i n
save data.mat