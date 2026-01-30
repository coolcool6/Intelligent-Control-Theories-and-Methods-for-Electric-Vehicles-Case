function [sys_pluse,N]=Generate_Pluse_Data(ny,nu)
%------获得脉冲响应系数-------%
% 采用单位阶跃响应的结果去算(书公式4.8)

sim_result=sim("get_step_data_2.slx");
time=sim_result.tout;
for i=1:length(sim_result.step11)
    if i==1
        pluse11(1)=sim_result.step11(1);
        pluse12(1)=sim_result.step12(1);
        pluse21(1)=sim_result.step21(1);
        pluse22(1)=sim_result.step22(1);
    else
        pluse11(i)=sim_result.step11(i)-sim_result.step11(i-1);
        pluse12(i)=sim_result.step12(i)-sim_result.step12(i-1);
        pluse21(i)=sim_result.step21(i)-sim_result.step21(i-1);
        pluse22(i)=sim_result.step22(i)-sim_result.step22(i-1);
    end
end

% ny=2;
% nu=2;
% %画图
% figure(1)
% subplot(2,2,1)
% plot(time,pluse11,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出1对输入1','FontSize',18);
% hold on
% subplot(2,2,2)
% plot(time,pluse12,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出1对输入2','FontSize',18);
% hold on
% subplot(2,2,3)
% plot(time,pluse21,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出2对输入1','FontSize',18);
% hold on
% subplot(2,2,4)
% plot(time,pluse22,"b",LineWidth=1);
% xlabel("时间/s",'FontSize',14);
% ylabel("响应",'FontSize',14);
% title('输出2对输入2','FontSize',18);


sys_pluse=zeros(length(pluse11),ny,nu);
sys_pluse(:,1,1)=pluse11;
sys_pluse(:,1,2)=pluse12;
sys_pluse(:,2,1)=pluse21;
sys_pluse(:,2,2)=pluse22;

size_pluse_data=size(sys_pluse);
sys_pluse=sys_pluse(2:size_pluse_data(1),:,:); %去掉0,截取剩余部分

N=25; %系统的过渡过程时间为N个采样间隔(输出图像后手动设置)
end
