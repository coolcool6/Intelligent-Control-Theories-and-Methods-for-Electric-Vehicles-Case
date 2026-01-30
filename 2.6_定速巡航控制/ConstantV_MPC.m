clear
clc
close all

Np = 5;
Ad = 1;
Bd = 0.1;

M_xx = zeros(Np, 1);
M_xu = zeros(Np, Np);
for i = 1:Np
    M_xx( i,: ) = Ad^i;
    for j = 1:Np
        if i >= j
            M_xu( i,j )= Ad^(i-j)*Bd;
        end
    end
end


x0 = [0];
N_max = 601;
x_sim = zeros(N_max+1,1);
x_sim(1,:) = x0';
u_sim = zeros(N_max,1);

x_min = 0;
x_max = 40;
u_min = -4;
u_max = 4;

u0 = repmat([0],[Np 1]);
Aeq = [];beq = [];
lb = [];ub = [];

for time = 1:300
    
    A = [ M_xu;
         -M_xu];              
    b = [ x_max - M_xx*x0;
         -x_min + M_xx*x0];

    options = optimoptions(@fmincon,'Display','Final','Algorithm','sqp',...
        'MaxIterations',300,'MaxFunEvals',30000);
    [U_best(:,time),fval(:,time),exitflag(:,time),output(:,time),...
        lambda(:,time),grad(:,time),hessian(:,:,time)] ...
                = fmincon('J_cost_cV35',u0,[],[],Aeq,beq,u_min*ones(Np,1),u_max*ones(Np,1),[],options);
    u_sim(time,:) = U_best(1,time)';
    x_sim(time+1,:) = (Ad*x0 + Bd*u_sim(time,:)')';
    x0 = x_sim(time+1,:)';

end

for time = 301:N_max
    
    A = [ M_xu;
         -M_xu];              
    b = [ x_max - M_xx*x0;
         -x_min + M_xx*x0];

    options = optimoptions(@fmincon,'Display','Final','Algorithm','sqp',...
        'MaxIterations',300,'MaxFunEvals',30000);
    [U_best(:,time),fval(:,time),exitflag(:,time),output(:,time),...
        lambda(:,time),grad(:,time),hessian(:,:,time)] ...
                = fmincon('J_cost_cV30',u0,[],[],Aeq,beq,u_min*ones(Np,1),u_max*ones(Np,1),[],options);
    u_sim(time,:) = U_best(1,time)';
    x_sim(time+1,:) = (Ad*x0 + Bd*u_sim(time,:)')';
    x0 = x_sim(time+1,:)';

end

figure(1)
subplot(1,2,1)
plot([0 0.1*N_max-1],[u_max u_max],'k--','LineWidth',0.5); hold on;
plot([0 0.1*N_max-1],[u_min u_min],'k--','LineWidth',0.5); hold on;
plot((0:0.1:0.1*(N_max-1)),u_sim(:,1),'r-','LineWidth',1.5);hold on;
xlabel('Time [s]')
ylabel('a [m/s^2]')
ylim([-4.5 4.5])
% xlim([0,15])


v_ref(1:301,1)=35;
v_ref(302:602,1)=15;


subplot(1,2,2)
% plot([0 0.1*N_max-1],[x_max x_max],'k--','LineWidth',0.5); hold on;
% plot([0 0.1*N_max-1],[x_min x_min],'k--','LineWidth',0.5); hold on;
plot((0:0.1:0.1*(N_max-1)),v_ref(1:end-1,1),'k-.','LineWidth',0.5,'DisplayName','v_{ref}');hold on;
plot((0:0.1:0.1*(N_max-1)),x_sim(1:end-1,1),'r-','LineWidth',1.5,'DisplayName','v');hold on;
xlabel('Time [s]')
ylabel('v [m/s]')
legend
ylim([0 40])
% xlim([0,15])


