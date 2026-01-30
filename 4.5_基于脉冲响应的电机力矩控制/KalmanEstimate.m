function [x_out,P_out,K]=KalmanEstimate(A,B,H,x_in,u_in,y_m,P_in,ny)
%本例中:
%状态量x=delta_Y_hat
%控制量u=delta_u
%测量量y_m=delta_ym
% A=Mhs,B=Hu,H=Cm*Ck

%获取维度
dix=length(x_in); %状态量的维度
diy_m=length(y_m); %获取观测量的维度

%模型误差和测量误差的协方差矩阵
Q_p=1*eye(dix);
R_p=0.01*eye(diy_m);

%卡尔曼滤波核心部分
x_pre=A*x_in+B*u_in;
P_pre=A*P_in*A'+Q_p;
K=P_pre*H'*inv(H*P_pre*H'+R_p);
x_out=x_pre+K*(y_m-H*x_pre);
P_out=(eye(dix)-K*H)*P_pre;
end