function [obj, x, rtime, gap, rtimeMP, rtimeSP, iter_num] = BDmaster_LDRsUB_AmbigSetProj(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                                      gamma1, gamma2, Mean, Cov, bigM, par)

% disp('------------------DRBP_Proj-----------------------------------');

[l, k] = size(W);

Omega = Cov + Mean*Mean';
[m, n] = size(A);

x = binvar(d,1);
lambda = sdpvar(1,1); 
vartheta = sdpvar(1,1); 
p0 = sdpvar(m,1);
Y = sdpvar(n,k,'full'); 
y0 = sdpvar(n,1);
rho = sdpvar(k,d,'full');
omega = sdpvar(m,d,'full');
T = sdpvar(m,l,'full'); 
Lambda = sdpvar(m,l,'full'); 
theta = sdpvar(d,1); 

constraints = [lambda>=0, T(:)>=0, Lambda(:)>=0, theta(:)>=0];  

constraints = [constraints, x + t <= s, sum(x) <= r]; 

% (37f-37g)
constraints = [constraints, A'*p0 + lambda*c0 == v0];
constraints = [constraints, A'*T*W + lambda*C == V];
constraints = [constraints, T*h + p0 >= 0];

% (37h) 
constraints = [constraints, A*Y + Lambda*W == lambda*B0];
tmp = 0;
for i = 1:d
    tmp = tmp + [b1i; b2i{i}; b3i]*theta(i);
end
constraints = [constraints, Lambda*h - A*y0 + tmp + lambda*b0 >= 0];

% (37i-37r)
for i = 1:d
   constraints = [constraints,...
       p0 - (1-x(i))*bigM*ones([m,1]) <= omega(:,i) <=  p0 + (1-x(i))*bigM*ones([m,1])]; %#ok<*CHAIN,*AGROW>
   constraints = [constraints, -x(i)*bigM*ones([m,1]) <= omega(:,i) <= x(i)*bigM*ones([m,1])];
   constraints = [constraints,...
       (T*W)'*[b1i; b2i{i}; b3i] - (1-x(i))*bigM*ones([k,1]) <= rho(:,i) <= (T*W)'*[b1i; b2i{i}; b3i] + (1-x(i))*bigM*ones([k,1])];
   constraints = [constraints, -x(i)*bigM*ones([k,1]) <= rho(:,i) <= x(i)*bigM*ones([k,1])];
   constraints = [constraints, theta(i) <= lambda];
   constraints = [constraints, lambda - bigM*(1-x(i)) <= theta(i) <= bigM*x(i)];
end

tmp2 = 0;
for i = 1:d
    tmp2 = tmp2 + [b1i; b2i{i};b3i]'*omega(:,i);
end

constraints = [constraints, vartheta >= trace(Omega*B0'*T*W) + trace(Omega*C'*Y) + Mean'*(sum(rho,2)+(T*W)'*b0+B0'*p0+C'*y0+Y'*c0) + ...
                tmp2 + b0'*p0 + c0'*y0];
        
obj = w'*x + vartheta;

options = sdpsettings('verbose',par,'solver', 'gurobi', 'cachesolvers', 1, 'savesolveroutput',1);
options.gurobi.IntFeasTol = 1e-9;

optimize(constraints, obj, options);
rtimeMP = ans.solvertime;     rtime = ans.solvertime;
rtimeSP = 0;    

sqrt_Cov = real(sqrtm(Cov));
iter_num = 0;
while 1
    iter_num = iter_num + 1;
    
    LB = value(obj);
    calG = value(0.5*((B0'*T*W + C'*Y) + (B0'*T*W + C'*Y)'));
    calg = value(sum(rho,2)+(T*W)'*b0+B0'*p0+C'*y0+Y'*c0);
    calA = sqrt_Cov*calG*sqrt_Cov;
    calb = 2*sqrt_Cov*calG*Mean + sqrt_Cov*calg;
    calc = Mean'*calG*Mean + calg'*Mean + value(tmp2 + b0'*p0 + c0'*y0);
    [subprob_obj, runtime, hat_mu, hat_Omega] = BDsubprob_LDRsUB_AmbigSetProj(calA, calb, calc, gamma1, gamma2);
    rtimeSP = rtimeSP + runtime; 
    rtime = rtime + runtime; 
    UB = w'*value(x) + subprob_obj;
    
    bar_mu = sqrt_Cov*hat_mu + Mean;
    bar_Omega = sqrt_Cov*hat_Omega*sqrt_Cov + Mean*hat_mu'*sqrt_Cov + (Mean*hat_mu'*sqrt_Cov)' + Mean*Mean';
    if abs(value(vartheta)-subprob_obj) <= 1e-6
        break;
    else
        constraints = [constraints, vartheta >= trace(bar_Omega*B0'*T*W) + trace(bar_Omega*C'*Y) + bar_mu'*(sum(rho,2)+(T*W)'*b0+B0'*p0+C'*y0+Y'*c0) + ...
                        tmp2 + b0'*p0 + c0'*y0];
    end
    
    optimize(constraints, obj, options);
    rtimeMP = rtimeMP + ans.solvertime;
    rtime = rtime +  ans.solvertime;
end    
    
gap = (UB-LB)/abs(LB)*100;
if round(UB) == 0 && round(LB) == 0
    gap = 0;  
end
x = value(x); obj = value(obj); %#ok<*NOANS>


