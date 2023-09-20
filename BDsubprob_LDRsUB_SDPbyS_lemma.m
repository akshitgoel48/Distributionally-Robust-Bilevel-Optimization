function [eta, sigmma, obj, status, time] = BDsubprob_LDRsUB_SDPbyS_lemma(x, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h, ...
                                                                          gamma1, gamma2, Mean, Cov, par)
% PRIMAL SUBPROBLEM      

[d, ~] = size(x);
[m, n] = size(A);
[l, k] = size(W);

varr = sdpvar(1,1);
vart = sdpvar(1,1);
Q = sdpvar(k,k); % symmetric  
M = sdpvar(k+1,k+1); % symmetric
q = sdpvar(k,1);
lambda = sdpvar(1,1);
tau = sdpvar(l,1);
p0 = sdpvar(m,1);
Y = sdpvar(n,k,'full');
y0 = sdpvar(n,1);
T = sdpvar(m,l,'full');
Lambda = sdpvar(m,l,'full');     

% non-negativity constraints
constraints = [lambda>=0, tau(:)>=0, T(:)>=0, Lambda(:)>=0]; 

% (42)   
tmp2 = 0; tmp3 = 0;
for i = 1:d
    tmp3 = tmp3 + (T*W)'*[b1i; b2i{i}; b3i]*x(i);
    tmp2 = tmp2 + [b1i; b2i{i};b3i]'*p0*x(i);
end
tmp1 = 0.5*(q - B0'*p0 - tmp3 - (T*W)'*b0 - C'*y0 - Y'*c0 - W'*tau);
constraints = [constraints,... 
                [Q  - 0.5*(B0'*T*W + C'*Y) - 0.5*(B0'*T*W + C'*Y)', tmp1; tmp1', varr - tmp2 - b0'*p0 - c0'*y0 + tau'*h] == M];
constraints = [constraints, (M >= 0):'psd2']; %#ok<*BDSCA>
%--------------------------------------------------------------------------

% (37c)
constraints = [constraints, vart >= sum(sum((gamma2*Cov + Mean*Mean').*Q)) + Mean'*q + ...
                sqrt(gamma1)*norm(sqrtm(Cov)*(q+2*Q*Mean))];
% (37d)
constraints = [constraints, (Q >= 0):'psd1'];

% (37f-37g)
constraints = [constraints, (A'*p0 + lambda*c0 == v0):'37f'];  
constraints = [constraints, (A'*(T*W) + lambda*C == V):'37g1'];
constraints = [constraints, (T*h + p0 >= 0):'37g2'];

% (41a) 
constraints = [constraints, (A*Y + Lambda*W == lambda*B0):'41a1'];
tmp = 0;
for i = 1:d    
    tmp = tmp + lambda*[b1i; b2i{i}; b3i]*x(i);
end
constraints = [constraints, (Lambda*h - A*y0 + tmp + lambda*b0 >= 0):'41a2'];

obj = varr + vart;        
options = sdpsettings('verbose',par,'solver','mosek','cachesolvers',1);

optimize(constraints, obj, options);
status = ans.info; %#ok<*NOANS>

if ~(contains(status, 'Successfully solved', 'IgnoreCase', true))
    options.mosek.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
    optimize(constraints, obj, options);   
end

obj = value(obj);
status = ans.info;
time = ans.solvertime;

DM = dual(constraints('psd2'));
eta = DM(1:k,k+1);
sigmma = dual(constraints('41a2'));
