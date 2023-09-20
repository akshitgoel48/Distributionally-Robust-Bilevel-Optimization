function [UB, coeff1, coeff2, coeff3, time, status] = BDsubprob_LDRsUB_IA_COP(x, d, w, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                                              gamma1, gamma2, Mean, Cov)
                                                  
% disp('------------------Subprob Inner Approx COP-----------------------------------');

[m, n] = size(A);
[l, k] = size(W);

varr = sdpvar(1,1);
vart = sdpvar(1,1);
q = sdpvar(k,1);
lambda = sdpvar(1,1); 
p0 = sdpvar(m,1);
Q = sdpvar(k,k);
Z = sdpvar(l,l,'full');
Y = sdpvar(n,k,'full'); 
y0 = sdpvar(n,1);
T = sdpvar(m,l,'full'); 
Lambda = sdpvar(m,l,'full'); 

% non-negativity constraints
Ztemp = 0.5*(Z + Z');
constraints = [lambda>=0, T(:)>=0, Lambda(:)>=0, Ztemp(:)>=0];    

% (37b,d)
tmp2 = 0;
tmp3 = 0;
for i = 1:d
    tmp2 = tmp2 + [b1i; b2i{i};b3i]'*p0*x(i);
    tmp3 = tmp3 +(T*W)'*[b1i; b2i{i}; b3i]*x(i);
end

tmp1 = 0.5*(q - B0'*p0 - tmp3 - (T*W)'*b0 - C'*y0 - Y'*c0);

constraints = [constraints, Q >= 0];

M = [Q  - 0.5*(B0'*T*W + C'*Y) - 0.5*(B0'*T*W + C'*Y)', tmp1; tmp1', varr - tmp2 - b0'*p0 - c0'*y0];    

% constraints = [constraints, (M == calH'*Ztemp*calH):'IA-COP'];

constraints = [constraints, (M(1:k,1:k) == W'*Ztemp*W)];
constraints = [constraints, (M(1:k,k+1) == -W'*Ztemp*h):'eta'];
constraints = [constraints, (M(k+1,k+1) == h'*Ztemp*h)];

% (37c)
constraints = [constraints,...
    vart >= sum(sum((gamma2*Cov + Mean*Mean').*Q)) + Mean'*q + sqrt(gamma1)*norm(sqrtm(Cov)*(q+2*Q*Mean))]; 

% (37f-37g)
constraints = [constraints, A'*p0 + lambda*c0 == v0 ];  % creating problem
constraints = [constraints, A'*(T*W) + lambda*C == V];
constraints = [constraints, T*h + p0 >= 0];

% (37h) 
constraints = [constraints, A*Y + Lambda*W == lambda*B0];
tmp = 0;
for i = 1:d
    tmp = tmp + lambda*[b1i; b2i{i}; b3i]*x(i);
end
constraints = [constraints, (Lambda*h - A*y0 + tmp + lambda*b0 >= 0):'sigma']; %#ok<*BDSCA>

obj = varr + vart;    
 
options = sdpsettings('verbose',0,'solver','mosek','cachesolvers',1,'savesolveroutput',1);
% 

solver_out = optimize(constraints, obj, options);

status = solver_out.info;

if ~(contains(status, 'Successfully solved', 'IgnoreCase', true))
    options.mosek.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF';
    solver_out = optimize(constraints, obj, options);
    status = solver_out.info;
end    

obj = value(obj);   time = solver_out.solvertime; 

sigmma = dual(constraints('sigma'));
eta = -dual(constraints('eta'))/2;

coeff1 = zeros(d,1);
coeff2 = zeros(d,1);
for i = 1:d
    bi = [b1i; b2i{i}; b3i];
    coeff1(i) = (-min(-value(sigmma)'*bi,0)+sum(max(bi,0)-min(bi,0))+...
                sum(max(value(eta),0)-min(value(eta),0)))*(1-x(i));
    coeff2(i) = (max(value(sigmma)'*bi,0)+sum(max(bi,0)-min(bi,0))+...
                sum(max(value(eta),0)-min(value(eta),0)))*x(i);    
end
coeff3 = obj;
UB = w'*x+obj;
