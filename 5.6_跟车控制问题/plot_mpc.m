function plot_mpc(u,xr,t)

figure(1)
plot(t,u,'LineWidth',1.5);
xlabel('时间（s）');
ylabel('u');
title('期望加速度输入');
hold on

figure(2)
plot(t,xr(1,1:end-1),'LineWidth',1.5);
xlabel('时间（s）');
ylabel('相对速度（m/s）');
title('两车相对速度');
hold on

figure(3)
plot(t,xr(2,1:end-1),'LineWidth',1.5);
xlabel('时间（s）');
ylabel('相对距离（m）');
title('两车相对距离');
hold on

figure(4)
plot(t,xr(3,1:end-1),'LineWidth',1.5);
xlabel('时间（s）');
ylabel('自车速度（m/s）');
title('自车速度');
hold on

end