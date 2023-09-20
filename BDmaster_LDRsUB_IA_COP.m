function [LB, xval, time, gap, master_time, subprob_time, numIters] = BDmaster_LDRsUB_IA_COP(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                                                             gamma1, gamma2, Mean, Cov, bigM)
                                     
x = binvar(d,1);
nu = sdpvar(1,1);
constraints = [sum(x) <= r, x + t <= s];

obj = w'*x + nu;
options = sdpsettings('verbose',0,'solver','gurobi','cachesolvers',1,'savesolveroutput',1);
options.gurobi.Presolve = 0;
options.gurobi.IntFeasTol = 1e-9;

xval = NaN*ones(d,1);   
LB = -inf;  UB = inf; eps = 1e-3;
time = 0;   flag = -1;  numIters = 0;
master_time = 0;    subprob_time = 0;

while ~(UB-LB < eps && flag == 0)
    
xval_previous = xval;

optimize(constraints, obj, options);

flag = 0;   solver_gap = NaN; %#ok<*NASGU>
if contains(ans.info, 'Unbounded', 'IgnoreCase', true)
    flag = -1;
elseif contains(options.solver, 'gurobi', 'IgnoreCase', true)
    solver_gap = ans.solveroutput.result.mipgap; %#ok<*NOANS>
end

xval = value(x); LB = value(obj); rtime = ans.solvertime; time = time + rtime;
master_time = master_time + rtime;

if round(xval_previous) == round(xval)
    % disp('Abnormal termination!!');
    break;
end

[UB_new, coeff1, coeff2, coeff3, rtime, ~] = BDsubprob_LDRsUB_IA_COP(xval, d, w, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                                          gamma1, gamma2, Mean, Cov);
                                     
if UB_new < UB
    UB = UB_new; end                     

constraints = [constraints, nu >= coeff3-bigM*coeff1'*x-bigM*coeff2'*(1-x)];                                      %#ok<*AGROW>

time = time + rtime;      subprob_time = subprob_time + rtime;                    

numIters = numIters + 1;

end

gap = (UB-LB)/abs(LB)*100;

if round(UB) == 0 && round(LB) == 0
    gap = 0;  
end