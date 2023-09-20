function [LB, x1, rtime, gap, master_time, subprob_time, iter_num] = BDmaster_LDRsUB_SDPbyS_lemma(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                                                                  gamma1, gamma2, Mean, Cov, bigM, eps1, time_limit, solver, par)

% disp('------------------BENDERS---------------------------------------');

warning('off','all')
        
x = binvar(d,1);
nu = sdpvar(1,1);
constraints = [sum(x) <= r, x + t <= s];

obj = w'*x + nu;
options = sdpsettings('verbose',0,'solver',solver,'cachesolvers',1,'savesolveroutput',1);
options.gurobi.Presolve = 0;
options.gurobi.IntFeasTol = 1e-9;

LB = -inf; UB = inf;
x1 = zeros(d,1);
rtime = 0; iter_num = 0;
master_time = 0;    subprob_time = 0;

while 1
   
    [eta, sigmma, obj1, status, time] = BDsubprob_LDRsUB_SDPbyS_lemma(x1, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h, ...
                                                                      gamma1, gamma2, Mean, Cov, par);     %#ok<ASGLU>
    rtime = rtime + time;
    subprob_time = subprob_time + time;
    iter_num = iter_num + 1;
 
    if w'*x1 + obj1 < UB 
        UB =  w'*x1 + obj1; end
    
    tmp1 = zeros(d,1);
    tmp2 = zeros(d,1);

    for i=1:d
        b_i = [b1i; b2i{i}; b3i];        
        tmp1(i) = x1(i)*(max(-b_i'*sigmma,0) + sum(max(b_i,0)-min(b_i,0))+ sum(max(eta,0)-min(eta,0)));
        tmp2(i) = (1-x1(i))*(-min(-b_i'*sigmma,0) + sum(max(b_i,0)-min(b_i,0)) + sum(max(eta,0)-min(eta,0)));  
    end
    
    diff = UB-LB;
    if diff > eps1
        al = obj1-bigM*sum(tmp1);
        ul = bigM*(tmp1-tmp2)';
        constraints = [constraints, nu >= ul*x+al]; %#ok<*AGROW> % optimality cut
    else              
        break;
    end
    
    % MASTER PROBLEM   
    optimize(constraints, obj, options);
    x1 = round(value(x)); LB = value(obj); rtime = rtime + ans.solvertime;
    master_time = master_time + ans.solvertime;
    solver_gap = ans.solveroutput.result.mipgap; %#ok<*NOANS,*NASGU>   
    
    if rtime > time_limit          
        break;
    end
end

gap = (UB-LB)/abs(LB)*100;

if round(UB) == 0 && round(LB) == 0
    gap = 0;  
end


