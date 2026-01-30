function J=J_cost_cV35(U)
    M_xx = evalin('base','M_xx');
    M_xu = evalin('base','M_xu');
    x0 = evalin('base','x0');
    Np = evalin('base','Np');
    X_k = M_xx * x0 + M_xu * U;
    J = sum(1*(X_k-35).^2 + 1*U.^2);
   
end
















