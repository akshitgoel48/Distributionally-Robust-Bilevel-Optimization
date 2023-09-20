function [obj, x, Y, y0, rtime] = Yanikoglu_Kuhn_PessTSBP(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                              Lb, xi, Mean, Cov, bigM, par, SLDR)

% disp('------------------TS-BP-----------------------------------');

[l, k] = size(W);

if SLDR == 1
    R = zeros(k, 2*k);    R(1,1) = 1;
    W = zeros(3*k,2*k);    h = zeros(3*k,1);
    xi = lift_space(xi, Mean);
    tempMean = mean(xi,2);       tempCov = cov(xi');  
    for i = 1:k
        R(i, 2*i-1:2*i) = 1;
        W(3*i-2, 2*i-1) = -1; h(3*i-2) = -Mean(i);
        W(3*i-1, 2*i) = 1; h(3*i-1) = 0;
        W(3*i, 2*i-1) = 1; W(3*i, 2*i) = -1; h(3*i) = Lb;                   
    end
    C = C*R;    V = V*R;    B0 = B0*R;
    Mean = tempMean; Cov = tempCov; 
    [l, k] = size(W);
end

Omega = Cov + Mean*Mean';
[m, n] = size(A);

x = binvar(d,1);
lambda = sdpvar(1,1); 
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
        
obj = w'*x + trace(Omega*B0'*T*W) + trace(Omega*C'*Y) + Mean'*(sum(rho,2)+(T*W)'*b0+B0'*p0+C'*y0+Y'*c0) + ...
        tmp2 + b0'*p0 + c0'*y0;

options = sdpsettings('verbose',par,'solver', 'gurobi', 'cachesolvers', 1, 'savesolveroutput',1);
options.gurobi.IntFeasTol = 1e-9;

optimize(constraints, obj, options);

x = value(x); obj = value(obj); rtime = ans.solvertime; %#ok<*NOANS>

Y = value(Y)/value(lambda);
y0 = value(y0)/value(lambda);

