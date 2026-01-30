%%主函数
% 状态转移矩阵考虑常数项的情况
% 阶跃响应的输出为id,iq,Te作为被控输出,
% 相比V3版本,测量输出增加id
% 即ny=2,nu=2,Cc=C
clear;
clc;
%%************* 1.设置参数和初始化 **************%%
%状态量:id,iq
%输入量:ud,uq
%控制输出:Te
%测量量:id,Te

%-------设置电机参数-------%
Ld=0.0085;
Lq=0.0085;
y_ref=2.875; %定子电阻
Pn=4; %电机极对数
Omega=100; %电机机械转速(单位:rad/s)
Fai=0.175; %转子磁链

%-------状态空间方程-------%
% 定义原始系统的状态空间模型,使用连续形式
A = [-y_ref/Ld Pn*Omega*Lq/Ld;-Pn*Omega*Ld/Lq -y_ref/Lq];
B = [1/Ld 0;0 1/Lq];
C = [0 3*Pn*Fai/2];
D = [0 0];
E = [0;-Pn*Omega*Fai/Lq]; %状态转移矩阵中的常数项
ts=0.001; %采样时间
sys_sta = ss(A, B, C, D);
Cc=C; %输出到控制输出的矩阵,和C一样
Cm=[1 0;0 3*Pn*Fai/2];

ny=2; %系统输出个数
nu=2; %控制输入个数
nx=2; %状态量个数
nyc=1; %系统控制输出个数
nym=2; %系统测量量个数

%-------DMC控制参数----------%
p=6; %预测时域
m=5; %控制时域

%跟踪参考值和控制输入的权重
w_y=diag(100*ones(1,p*nyc));
w_u=diag(0.1*ones(1,m*nu));

% w_y_diag=zeros(1,p*nyc);
% for i=1:nyc:(p*nyc)
%     if i==1
%         w_y_diag(1,1:nyc)=100;
%     else
%         w_y_diag(1,i:(i+nyc-1))=w_y_diag(i-nyc)*0.9;
%     end
% end
% w_y=diag(w_y_diag);
% 
% w_u_diag=zeros(1,m*nu);
% for i=1:nu:(m*nu)
%     if i==1
%         w_u_diag(1,1:nu)=10;
%     else
%         w_u_diag(1,i:(i+nu-1))=w_u_diag(i-nu)*0.9;
%     end
% end
% w_u=diag(w_u_diag);

%------得到阶跃响应的结果-------%
[sys_step,N]=Generate_Step_Data(ny,nu);

%-------计算Kmpc和相关矩阵---------%
[Kmpc,Mss,Su,M_hua,Su_hua,Ck]=DMC_Matrix_Calculate(sys_step,N,ny,nu,p,m,Cc,w_y,w_u);

%-------设置仿真参数--------%
T_sim=2; %仿真时长
k_sim=ts:ts:T_sim; %仿真步
step_total=length(k_sim); %总仿真步数
delta_u=zeros(nu,step_total);
u=zeros(nu,step_total); %输入值
ym=zeros(nym,step_total); %测量值
y_ref=zeros(1,step_total); %参考值
Y_hat=zeros(N*ny,step_total); %Y估计值
Ep=zeros(p,step_total); %误差
R=zeros(p,step_total); %参考值序列

% P_kalman=zeros(ny*N,ny*N); %卡尔曼估计误差协方差矩阵初值
P_kalman=0.1*ones(ny*N,ny*N);
H_kalman=Cm*Ck; %卡尔曼滤波估计的测量矩阵


x=zeros(nx,step_total); %状态量
y=zeros(1,step_total); %控制输出量


% % 设置仿真初始值2,正弦信号的
% delta_u0=zeros(nu,1);
% u0=zeros(nu,1);
% ym0=0;
% Y_hat0=zeros(N*ny,1);
% 
% %设置参考值2，正弦信号
% f = 10; % 频率
% Amp = 20; % 振幅
% phi = 0; % 相位
% y_ref=Amp * sin(2*pi*f*k_sim + phi);

% % 设置仿真初始值3,连续变化信号的
% delta_u0=zeros(nu,1);
% u0=zeros(nu,1);
% ym0=zeros(nym,1);
% Y_hat0=zeros(N*ny,1);
% 
% %设置参考值3,连续变化信号
% y_ref(1,1500:3500)=10;
% y_ref(1,3501:5000)=10:0.01:25-0.01;
% y_ref(1,5001:7000)=25;
% y_ref(1,7001:9000)=15;
% y_ref(1,9001:11000)=5;
% y_ref(1,11001:13000)=20;
% y_ref(1,13001:step_total)=10;

% 设置仿真初始值4,连续变化信号的
delta_u0=zeros(nu,1);
u0=zeros(nu,1);
ym0=zeros(nym,1);
Y_hat0=zeros(N*ny,1);

%设置参考值4,连续变化信号
y_ref(1,150:350)=10;
y_ref(1,351:500)=10:0.1:25-0.1;
y_ref(1,501:700)=25;
y_ref(1,701:900)=15;
y_ref(1,901:1100)=5;
y_ref(1,1101:1300)=20;
y_ref(1,1301:step_total)=10;


Ad=A*ts+eye(nx);
%-----离散状态空间模型-----%
Bd=B*ts;
Cd=C;
Dd=D;
Ed=E*ts;

%%********进入仿真大循环**********%%
for k=1:step_total

    %%***** 2.状态估计(使用卡尔曼滤波) *******%%
    ym(:,k)=Cm*x(:,k); %获取当前测量值,注意状态量对应的是DMC里的系统输出
    if k==1
        [Y_hat(:,k),P_temp,K]=KalmanEstimate(Mss,Su,H_kalman,Y_hat0,delta_u0,ym(:,k),P_kalman,ny);
    else
        [Y_hat(:,k),P_temp,K]=KalmanEstimate(Mss,Su,H_kalman,Y_hat(:,k-1),delta_u(:,k-1),ym(:,k),P_kalman,ny);
    end
    P_kalman=P_temp; %更新协方差矩阵

    %%***** 3.误差计算和得到该时刻控制量 ******%% 
    %参考值序列
    for i=1:p
        if ((k+i)<=step_total)
            R(i,k)=y_ref(k+i);
        else
            R(i,k)=y_ref(step_total);
        end
    end
    Ep(:,k)=R(:,k)-M_hua*Y_hat(:,k);  %误差计算
    
    %%****** 4.该时刻控制量计算 ********%%
    delta_u(:,k)=Kmpc*Ep(:,k);
    if k==1
        u(:,k)=u0+Kmpc*Ep(:,k);
    else
        u(:,k)=u(:,k-1)+Kmpc*Ep(:,k);
    end

    %%***** 5.得到该时刻实际输出值和测量值 *****%%
    x(:,k+1)=Ad*x(:,k)+Bd*u(:,k)+Ed; %下一时刻状态量
    y(:,k)=Cd*x(:,k)+Dd*u(:,k);
end

%%绘图程序
figure(1)
plot(k_sim,y(1,:),"b",LineWidth=1);
hold on
plot(k_sim,y_ref(1,:),"r--",LineWidth=1);
xlabel("时间/s",'FontSize',14);
ylabel("转矩/Nm",'FontSize',14);
legend('实际转矩','目标转矩','FontSize',12);
title('转矩跟踪效果图像','FontSize',18);