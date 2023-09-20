function [LB, xval, time, gap, master_time, subprob_time, numIters] = BDmaster_DiscrtLB_SDP(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W,... 
                                                                                            gamma1, gamma2, xi0, xi, bigM, par, presolve_subprob)

% disp('------------------Lower bound Discrete BD------------------');

[m, ~] = size(A);

nsamp = size(xi,2);

%------------------------------------------------------------------
x = binvar(d,1);
nu = sdpvar(1,1); 

%------------------------------------------------------------------
obj = w'*x+nu; 
constraints = [x + t <= s, sum(x) <= r]; 

options = sdpsettings('verbose',par,'solver','gurobi','cachesolvers',1,'savesolveroutput',1);
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
if contains(ans.info, 'Unbounded', 'IgnoreCase', true) %#ok<*NOANS>
    flag = -1;
elseif contains(options.solver, 'gurobi', 'IgnoreCase', true)
    solver_gap = ans.solveroutput.result.mipgap;
end

xval = value(x); LB = value(obj); rtime = ans.solvertime; time = time + rtime;
master_time = master_time + rtime;

if round(xval_previous) == round(xval)
    % disp('Abnormal termination!!');
    break;
end

[UB_new, cut, gamma, sigma, rtime, ~] = BDsubprob_DiscrtLB_SDP(xval, d, w, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W,... 
                                                               gamma1, gamma2, xi0, xi, 0, presolve_subprob);
if UB_new < UB
    UB = UB_new; end    

time = time + rtime;             
subprob_time = subprob_time + rtime;

numIters = numIters + 1;

tmp = zeros(d,1);
for i = 1:d
    for s = 1:nsamp
        tmp(i) = tmp(i)-sigma(:,s)'*[b1i; b2i{i}; b3i];
    end
end

xval = round(xval); 
for i= 1:d
    if xval(i)==1
        cut = cut - bigM*(1-x(i))*max(tmp(i),0);        
        for s = 1:nsamp
            tmp2 = gamma(s)*[b1i; b2i{i}; b3i];
            cut = cut - bigM*(1-x(i))*(max(-tmp2,0)+max(tmp2,0))'*ones(m,1);
        end
    else
       cut = cut - bigM*x(i)*max(-tmp(i),0);        
        for s = 1:nsamp
            tmp2 = gamma(s)*[b1i; b2i{i}; b3i];
            cut = cut - bigM*x(i)*(max(-tmp2,0)+max(tmp2,0))'*ones(m,1);
        end
    end
end

constraints = [constraints, nu >= cut]; %#ok<*AGROW>

end

gap = (UB-LB)/abs(LB)*100;

if round(UB) == 0 && round(LB) == 0
    gap = 0;  
end

