%% 预测矩阵生成，解决优化问题
function [P,S,W]=predmodgen(LTI,dim)

LTI.C=eye(dim.nx);
dim.ny=dim.nx;

%初始状态的预测矩阵
P=zeros(dim.ny*(dim.N),dim.nx);
for k=0:dim.N-1
    P(k*dim.ny+1:(k+1)*dim.ny,:)=LTI.C*LTI.A^k;
end

%输出的预测矩阵
S=zeros(dim.ny*(dim.N),dim.nu*(dim.N));
for k=1:dim.N-1
    for i=0:k-1
        S(k*dim.ny+1:(k+1)*dim.ny,i*dim.nu+1:(i+1)*dim.nu) = LTI.C*LTI.A^(k-1-i)*LTI.B;
    end
end

W = zeros(dim.ny,dim.nu*dim.N);
for i=0:dim.N-1
    W(1:dim.ny,i*dim.nu+1:(i+1)*dim.nu) = LTI.C*LTI.A^(dim.N-i-1)*LTI.B; %(i+1)*ones(dim.ny,dim.nu);
%     W(k*dim.ny+1:(k+1)*dim.ny,i*dim.nu+1:(i+1)*dim.nu) = LTI.C*LTI.A^(k-1-i)*LTI.B;
end