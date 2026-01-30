function [Kmpc,Mhs,Hu,M_hua,Hu_hua,Ck,I_hua]=MAC_Matrix_Calculate(sys_pluse,N,ny,nu,nc,p,m,Cc,w_y,w_u)
%-----计算不变输入非零初始状态响应的矩阵Mhs----%
Mhs_cell=cell(N,N);
for i=1:N
    for j=1:N
        if i<N
            if j==i+1
                Mhs_cell{i,j}=diag(ones(1,ny));
            else
                Mhs_cell{i,j}=zeros(ny,ny);
            end
        else
            Mhs_cell{i,j}=zeros(ny,ny);
        end
    end
end
Mhs=cell2mat(Mhs_cell);

%----计算零初始条件下,系统对输入的响应矩阵Hu----%
Hu_cell=cell(N,1);
for count=1:N
    temp=zeros(ny,nu);
    for i=1:nu
        for j=1:ny
            temp(j,i)=sys_pluse(count,j,i);
        end
    end
    Hu_cell{count}=temp;
end
Hu=cell2mat(Hu_cell);

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
%这里计算不同于DMC,计算逻辑还要修改下
M_hua_cell=cell(p,N);
for i=1:p
    for j=1:N
        if (j>=2) && (j<=i+1)
            M_hua_cell{i,j}=Cc;
        else
            M_hua_cell{i,j}=zeros(size(Cc));
        end
    end
end
M_hua=cell2mat(M_hua_cell);

%----计算花写Hu------%
Hu_hua_cell=cell(p,m);
size1=size(Cc*Hu_cell{1});
for i=1:p
    for j=1:m
        index=i-j+1;
        if index>0
            Hu_cell_sum=zeros(ny,nu);
            for count=1:1:index
                Hu_cell_sum=Hu_cell_sum+Hu_cell{index};
            end
            Hu_hua_cell{i,j}=Cc*Hu_cell_sum;
        else
            Hu_hua_cell{i,j}=zeros(size1);
        end
    end
end
Hu_hua=cell2mat(Hu_hua_cell);

%-----计算花写I-----%
I_hua_cell=cell(p,1);
for i=1:p
    I_hua_cell{i,1}=eye(nc);
end
I_hua=cell2mat(I_hua_cell);

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
Kmpc=K1*inv(Hu_hua'*(w_y'*w_y)*Hu_hua+w_u'*w_u)*Hu_hua'*(w_y'*w_y);
end