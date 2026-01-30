clear
clc
close all

Np = 5;
Ad = [1 -0.1; 1 -0.5];
Bd = [0.5 0; 0 -0.5];

M_xx = zeros(2*Np, 2);
M_xu = zeros(2*Np, 2*Np);
for i = 1:Np
    M_xx( (2*i-1:2*i),: ) = Ad^i;
    for j = 1:Np
        if i >= j
            M_xu( (2*i-1:2*i),(2*j-1:2*j) )= Ad^(i-j)*Bd;
        end
    end
end

T1=zeros(Np,2*Np);
T2=zeros(Np,2*Np);
for i=1:1:Np
    T1(i,2*i-1)=1;
    T2(i,2*i  )=1;
end

x0 = [3; -2];
N_max = 31;
x_sim = zeros(N_max+1,2);
x_sim(1,:) = x0';
u_sim = zeros(N_max,2);

x1_min = 2;
x2_min = -3;
x1_max = 9;
x2_max = 6;
u1_min = -2;
u2_min = -3;
u1_max = 2;
u2_max = 3;

u0 = repmat([0;0],[Np 1]);
Aeq = [];beq = [];
lb = [];ub = [];

for time = 1:N_max
    
    A = [ T1
         -T1
          T2
         -T2
          T1*M_xu;
         -T1*M_xu;
          T2*M_xu;
         -T2*M_xu];              
    b = [ u1_max*ones(Np,1)
         -u1_min*ones(Np,1)
          u2_max*ones(Np,1)
         -u2_min*ones(Np,1)
          x1_max - T1*M_xx*x0;
         -x1_min + T1*M_xx*x0;
          x2_max - T2*M_xx*x0;
         -x2_min + T2*M_xx*x0];

    options = optimoptions(@fmincon,'Display','Final','Algorithm','sqp',...
        'MaxIterations',300,'MaxFunEvals',30000);
    [U_best(:,time),fval(:,time),exitflag(:,time),output(:,time),...
        lambda(:,time),grad(:,time),hessian(:,:,time)] ...
                = fmincon('J_cost',u0,A,b,Aeq,beq,lb,ub,[],options);
    u_sim(time,:) = U_best(1:2,time)';
    x_sim(time+1,:) = (Ad*x0 + Bd*u_sim(time,:)')';
    x0 = x_sim(time+1,:)';

end

figure(1)
subplot(2,2,1)
plot([0 N_max-1],[u1_max u1_max],'k--','LineWidth',0.5); hold on;
plot([0 N_max-1],[u1_min u1_min],'k--','LineWidth',0.5); hold on;
plot((0:N_max-1),u_sim(:,1),'r-','LineWidth',1.5);hold on;
xlabel('Time [s]')
ylabel('u_1')
ylim([-2.5 2.5])

subplot(2,2,2)
plot([0 N_max-1],[u2_max u2_max],'k--','LineWidth',0.5); hold on;
plot([0 N_max-1],[u2_min u2_min],'k--','LineWidth',0.5); hold on;
plot((0:N_max-1),u_sim(:,2),'r-','LineWidth',1.5);hold on;
xlabel('Time [s]')
ylabel('u_2')
ylim([-4 4])

subplot(2,2,3)
plot([0 N_max-1],[x1_max x1_max],'k--','LineWidth',0.5); hold on;
plot([0 N_max-1],[x1_min x1_min],'k--','LineWidth',0.5); hold on;
plot((0:N_max-1),x_sim(1:end-1,1),'r-','LineWidth',1.5);hold on;
xlabel('Time [s]')
ylabel('x_1')
ylim([1 10])

subplot(2,2,4)
plot([0 N_max-1],[x2_max x2_max],'k--','LineWidth',0.5); hold on;
plot([0 N_max-1],[x2_min x2_min],'k--','LineWidth',0.5); hold on;
plot((0:N_max-1),x_sim(1:end-1,2),'r-','LineWidth',1.5);hold on;
xlabel('Time [s]')
ylabel('x_2')
ylim([-4 7])


