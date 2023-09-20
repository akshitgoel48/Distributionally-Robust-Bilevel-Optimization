function [obj_val, cut_init, gamma_dual, sigma_dual, rtime, status] = BDsubprob_DiscrtLB_SDP(x, d, w, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W,... 
                                                                                             gamma1, gamma2, xi0, xi, par, presolve)

x = round(x);

[m, n] = size(A);
[~, k] = size(W);

Mean = mean(xi0,2); % Mean0;
Cov = cov(xi0'); % Cov0;

nsamp = size(xi,2);

% PRIMAL
%------------------------------------------------------------------
lambda = sdpvar(1,1); 
p = sdpvar(m,nsamp,'full');
y = sdpvar(n,nsamp,'full');
varr = sdpvar(1,1);
vart = sdpvar(1,1);
Q = sdpvar(k,k); % symmetric  
q = sdpvar(k,1);

%------------------------------------------------------------------
obj = w'*x; 
constraints = []; 
% obj1 = 0;
for ns = 1:nsamp 
    c = C*xi(:,ns) + c0;
    v = V*xi(:,ns) + v0;    
    
    tmp2 = 0;
    for i = 1:d
        tmp2 = tmp2 + [b1i; b2i{i};b3i]'*p(:,ns)*x(i);
    end
    constraints = [constraints, (varr >= (B0*xi(:,ns))'*p(:,ns) + tmp2 + b0'*p(:,ns) + c'*y(:,ns)-xi(:,ns)'*Q*xi(:,ns)-xi(:,ns)'*q):strcat('gamma',ns)]; %#ok<*AGROW>
    constraints = [constraints, A'*p(:,ns) + lambda*c == v];
    tmp = 0;
    for i = 1:d
        tmp = tmp + [b1i; b2i{i}; b3i]*lambda*x(i);
    end
    constraints = [constraints, (A*y(:,ns) <= tmp + lambda*(B0*xi(:,ns)+b0)):strcat('sigma',ns)];
end

constraints = [constraints, vart >= sum(sum((gamma2*Cov + Mean*Mean').*Q)) + Mean'*q + sqrt(gamma1)*norm(sqrtm(Cov)*(q+2*Q*Mean))];
obj1 = varr + vart;

constraints = [constraints, [lambda>=0, p(:)>=0]]; % theta(:)>=0, 

constraints = [constraints, (Q >= 0):'psd1']; %#ok<*BDSCA>

options = sdpsettings('verbose',par,'solver','mosek');
if presolve == 0
    options.mosek.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
end

optimize(constraints, obj1, options);
%------------------------------------------------------------------
        
obj_val = value(obj+obj1); rtime = ans.solvertime; status = ans.info; %#ok<*NOANS>

gamma_dual = zeros(nsamp,1);
sigma_dual = zeros(m,nsamp);
for ns = 1:nsamp 
    gamma_dual(ns) = dual(constraints(strcat('gamma',ns)));
    sigma_dual(:,ns) = dual(constraints(strcat('sigma',ns)));
end
cut_init = value(obj1);
