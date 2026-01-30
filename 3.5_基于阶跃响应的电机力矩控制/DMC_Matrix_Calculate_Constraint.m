function [Mss,Su,M_hua,Su_hua,Ck,Mb_hua,Sub_hua,H_qp,A_qp]=DMC_Matrix_Calculate_Constraint(sys_step,N,ny,nu,p,m,Cc,w_y,w_u,Cb)
%-----计算不变输入非零初始状态响应的矩阵Mss----%
Mss_cell=cell(N,N);
for i=1:N
    for j=1:N
        if i<N
            if j==i+1
                Mss_cell{i,j}=diag(ones(1,ny));
            else
                Mss_cell{i,j}=zeros(ny,ny);
            end
        else
            if j==i
                Mss_cell{i,j}=diag(ones(1,ny));
            else
                Mss_cell{i,j}=zeros(ny,ny);
            end
        end
    end
end
Mss=cell2mat(Mss_cell);

%----计算零初始条件下,系统对输入的响应矩阵Su----%
Su_cell=cell(N,1);
for count=1:N
    temp=zeros(ny,nu);
    for i=1:nu
        for j=1:ny
            temp(j,i)=sys_step(count,j,i);
        end
    end
    Su_cell{count}=temp;
end
Su=cell2mat(Su_cell);

%----计算由定义响应Y到k时刻输出y的矩阵C----%
%----为和状态空间方程的C区分，这里用Ck-----%
Ck_cell=cell(1,N);
for i=1:N
    if i==1
        Ck_cell{i}=diag(ones(1,ny));
    else
        Ck_cell{i}=zeros(ny,ny);
    end
end
Ck=cell2mat(Ck_cell);

%-----计算花写M--------%
M_hua_cell=cell(p,N);
for i=1:p
    for j=1:N
        if j==i+1
            M_hua_cell{i,j}=Cc;
        else
            M_hua_cell{i,j}=zeros(size(Cc));
        end
    end
end
M_hua=cell2mat(M_hua_cell);

%----计算花写Su------%
Su_hua_cell=cell(p,m);
size1=size(Cc*Su_cell{1});
for i=1:p
    for j=1:m
        index=i-j+1;
        if index>0
            Su_hua_cell{i,j}=Cc*Su_cell{index};
        else
            Su_hua_cell{i,j}=zeros(size1);
        end
    end
end
Su_hua=cell2mat(Su_hua_cell);

%-----计算H---------%
H=Su_hua'*(w_y'*w_y)*Su_hua+w_u'*w_u;

%-----计算控制增量约束转换的T------%
T_cell=cell(m,m);
for i=1:m
    for j=1:m
        if i==j
            T_cell{i,j}=eye(nu);
        else
            T_cell{i,j}=zeros(nu,nu);
        end
    end
end
T=cell2mat(T_cell);

%-----计算控制量约束转换的L------%
L_cell=cell(m,m);
for i=1:m
    for j=1:m
        if j<=i
           L_cell{i,j}=eye(nu); 
        else
           L_cell{i,j}=zeros(nu,nu);
        end
    end
end
L=cell2mat(L_cell);

%-----计算输出约束转换矩阵Mb_hua和Sub_hua----%
%套路和算M_hua和Su_hua一样，Cc换成Cb
Mb_hua_cell=cell(p,N);
for i=1:p
    for j=1:N
        if j==i+1
            Mb_hua_cell{i,j}=Cb;
        else
            Mb_hua_cell{i,j}=zeros(size(Cb));
        end
    end
end
Mb_hua=cell2mat(Mb_hua_cell);

Sub_hua_cell=cell(p,m);
size2=size(Cb*Su_cell{1});
for i=1:p
    for j=1:m
        index=i-j+1;
        if index>0
            Sub_hua_cell{i,j}=Cb*Su_cell{index};
        else
            Sub_hua_cell{i,j}=zeros(size2);
        end
    end
end
Sub_hua=cell2mat(Sub_hua_cell);

%---- qp问题所需的离线矩阵------%
H_qp=2*H;
Cu=[-T',T',-L',L',-Sub_hua',Sub_hua']';
A_qp=-Cu;

end