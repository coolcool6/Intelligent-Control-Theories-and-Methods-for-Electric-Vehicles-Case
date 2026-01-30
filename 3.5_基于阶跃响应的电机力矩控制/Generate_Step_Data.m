function [sys_step,N]=Generate_Step_Data(ny,nu)
%------阶跃响应测试-------%

sim_result=sim("get_step_data.slx");
time=sim_result.tout;
step11=sim_result.step11;
step12=sim_result.step12;
step21=sim_result.step21;
step22=sim_result.step22;

% ny=2;
% nu=2;
% %画图
% figure(1)
% subplot(2,2,1)
% plot(time,step11,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出1对输入1','FontSize',18);
% hold on
% subplot(2,2,2)
% plot(time,step12,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出1对输入2','FontSize',18);
% hold on
% subplot(2,2,3)
% plot(time,step21,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出2对输入1','FontSize',18);
% hold on
% subplot(2,2,4)
% plot(time,step22,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出2对输入2','FontSize',18);


sys_step=zeros(length(step11),ny,nu);
sys_step(:,1,1)=step11;
sys_step(:,1,2)=step12;
sys_step(:,2,1)=step21;
sys_step(:,2,2)=step22;

size_step_data=size(sys_step);
sys_step=sys_step(2:size_step_data(1),:,:); %去掉0,截取剩余部分

N=25; %系统的过渡过程时间为N个采样间隔(输出图像后手动设置)
end
