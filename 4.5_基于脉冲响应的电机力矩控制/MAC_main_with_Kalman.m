%%主函数
% 状态转移矩阵考虑常数项的情况
% 修改阶跃响应的输出,为id,iq;Te作为被控输出,测量输出id,Te
% 即ny=2,nu=2,Cc=C,Cm=Cc
% 20240318:使用单位阶跃响应的系数做差,可以出来好结果
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
% Cm=[1 0;0 3*Pn*Fai/2];
Cm=[0 3*Pn*Fai/2];

ny=2; %系统输出个数
nu=2; %控制输入个数
nx=2; %状态量个数
nyc=1; %系统控制输出个数
nym=1; %系统测量量个数

%-------MAC控制参数----------%
p=8; %预测时域
m=7; %控制时域

%跟踪参考值和控制输入的权重(重要设置!!)
w_y=diag(1*ones(1,p*nyc));
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

%------得到脉冲响应的结果-------%
[sys_pluse,N]=Generate_Pluse_Data(ny,nu);

%-------计算Kmpc和相关矩阵---------%
[Kmpc,Mhs,Hu,M_hua,Hu_hua,Ck,I_hua]=MAC_Matrix_Calculate(sys_pluse,N,ny,nu,nyc,p,m,Cc,w_y,w_u);

%-------设置仿真参数--------%
T_sim=2; %仿真时长
k_sim=ts:ts:T_sim; %仿真步
step_total=length(k_sim); %总仿真步数
delta_u=zeros(nu,step_total);
u=zeros(nu,step_total); %输入值
ym=zeros(nym,step_total); %测量值
delta_ym=zeros(nym,step_total); %测量值变化量，用于状态估计
y_ref=zeros(1,step_total); %参考值
Y_hat=zeros(N*ny,step_total); %Y估计值
delta_Y_hat=zeros(N*ny,step_total); %Y估计值变化量，状态估计的结果，用于预测方程计算
Ep=zeros(p,step_total); %误差
R=zeros(p,step_total); %参考值序列
yc_hat=zeros(nyc,step_total); %被控输出估计值

P_kalman=0.1*ones(ny*N,ny*N);
H_kalman=Cm*Ck; %卡尔曼滤波估计的测量矩阵

% % P_kalman=zeros(ny*N,ny*N); %卡尔曼估计误差协方差矩阵初值
% P_kalman=0.1*ones(ny*N,ny*N);
% H_kalman=Cm*Ck; %卡尔曼滤波估计的测量矩阵

x=zeros(nx,step_total); %状态量
y=zeros(1,step_total); %控制输出量

% 设置仿真初始值4,连续变化信号的
delta_u0=zeros(nu,1);
u0=zeros(nu,1);
delta_ym0=zeros(nym,1);
ym0=zeros(nym,1);
delta_Y_hat0=zeros(N*ny,1);
Y_hat0=zeros(N*ny,1);
yc_hat0=zeros(nyc,1);

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

    %%***** 2.状态估计 *******%%
    ym(:,k)=Cm*x(:,k); %获取当前测量值,注意状态量对应的是MAC里的系统输出
    if k==1
        [delta_Y_hat(:,k),P_temp,K]=KalmanEstimate(Mhs,Hu,H_kalman,delta_Y_hat0,delta_u0,delta_ym0,P_kalman,ny);
        delta_ym(:,k)=ym(:,k)-ym0;
        yc_hat(:,k)=Cc*Ck*delta_Y_hat(:,k)+yc_hat0;
    else
        [delta_Y_hat(:,k),P_temp,K]=KalmanEstimate(Mhs,Hu,H_kalman,delta_Y_hat(:,k-1),delta_u(:,k-1),delta_ym(:,k-1),P_kalman,ny);
        delta_ym(:,k)=ym(:,k)-ym(:,k-1);
        yc_hat(:,k)=Cc*Ck*delta_Y_hat(:,k)+yc_hat(:,k-1);
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
    Ep(:,k)=R(:,k)-M_hua*delta_Y_hat(:,k)-I_hua*yc_hat(:,k);  %误差计算
    
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
ylabel("转矩/nm",'FontSize',14);
legend('实际转矩','目标转矩','FontSize',12);
title('转矩跟踪效果图像','FontSize',18);