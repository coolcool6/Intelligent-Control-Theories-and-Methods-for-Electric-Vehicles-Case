function [Kmpc,Mss,Su,M_hua,Su_hua,Ck]=DMC_Matrix_Calculate(sys_step,N,ny,nu,p,m,Cc,w_y,w_u)
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

%-----计算Kmpc------%
%跟踪参考值和控制输入的权重

% Kmpc表达式的第一部分
K1_cell=cell(1,m);
for i=1:m
    if i==1
        K1_cell{1,i}=diag(ones(1,nu));
    else
        K1_cell{1,i}=zeros(nu,nu);
    end
end
K1=cell2mat(K1_cell);
Kmpc=K1*inv(Su_hua'*(w_y'*w_y)*Su_hua+w_u'*w_u)*Su_hua'*(w_y'*w_y);
end