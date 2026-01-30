%% ACC-MPC结果
clc
clear
%% 系统方程
Teng=0.46;
Keng=0.732;
Af=-1/Teng;
Bf=-Keng/Teng;
Cf=eye(3);
Thw=1.6;
Ts=0.05;
T_total=10;
T=T_total/Ts;
%初始速度
v0=15;
%初始间距
init_dist=5;

%% 连续系统离散化
At=[0,1,-Thw;0,0,-1;0,0,Af];
Bt=[0;0;Bf];
plant=ss(At,Bt,Cf,0);
sys1=c2d(plant,Ts);
A=sys1.A;
B=sys1.B;
C=sys1.C;

%% 自车参数
%自车速度
vh=host_velocity(v0,T);
sec_row=[init_dist*ones(1,length(vh))];
x0=[sec_row;sec_row;vh];
u=zeros(1,length(vh));

xr=[zeros(1,length(vh));zeros(1,length(vh));zeros(1,length(vh))];

x=[zeros(1,T_total);zeros(1,T_total);zeros(1,T_total)];
s_lb=[0;0;15];
s_ub=[2;2.5;40];
%定义LTI系统
LTI.A=A;
LTI.B=B;
LTI.C=eye(3);

%定义dim系统
dim.nx=length(A);
dim.ny=length(B);
dim.nu=1;
dim.N=20;
N=dim.N;

%定义cost函数
Q=C'*C;
R=1;

%% 预测模型和cost函数
[P,S]=predmodgen(LTI,dim);
[H,h]=costgen(P,S,Q,R,dim);

%% MPC仿真
umin=-3*ones(1,N);
umax=5*ones(1,N);
xr(:,1)=x0(:,1);
y(:,1)=C*x0(:,1);
for i=1:T
    t(i)=(i-1)*Ts;
    f=h*xr(:,i);
    Aueq=C*x0(:,1:N);
    bueq=y(:,1);
    options=optimoptions('quadprog','Display','off');
    warning off;
    [ures,~,exitflag]=quadprog(H,f,[],[],Aueq,bueq,umin,umax,[],options);
    u(i)=ures(1);
    xr(:,i+1)=A*xr(:,i)+B*u(i);
    y(:,i)=C*xr(:,i+1);
    x(:,1)=xr(:,i);
        for j=1:N-1
            x(:,j+1)=A*x(:,j)+B*ures(j);
            y_res(:,j)=C*x(:,j);
        end
end

%% 绘制结果
plot_mpc(u,xr,t);