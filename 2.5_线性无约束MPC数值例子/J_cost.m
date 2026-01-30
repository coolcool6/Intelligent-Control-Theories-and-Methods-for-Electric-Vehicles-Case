function J=J_cost(U)
    M_xx = evalin('base','M_xx');
    M_xu = evalin('base','M_xu');
    x0 = evalin('base','x0');
    Np = evalin('base','Np');
    X_k = M_xx * x0 + M_xu * U;
    J = sum(X_k.^2-2*(repmat([8;7],[Np 1])).*X_k) + 0.1*sum(U.^2);
end
















