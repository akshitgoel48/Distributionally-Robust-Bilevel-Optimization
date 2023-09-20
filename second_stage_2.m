function obj =  second_stage_2(x, xi, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0)

m = size(A,1);
n = size(A,2);
d = size(x,1);

tmp1 = 0;         
for i = 1:d
    tmp1 = tmp1 + [b1i; b2i{i}; b3i]*x(i);
end
tmp1 = tmp1 + b0;

options = sdpsettings('verbose',0,'solver','gurobi','cachesolvers',1);

nsamp = size(xi,2);
    
lambda = sdpvar(1,nsamp); % nonnegative;
p = sdpvar(m,nsamp,'full'); % nonnegative;
y = sdpvar(n,nsamp,'full');

c = C*xi+c0;
v = V*xi+v0;  
       
tmp = tmp1 + B0*xi;

obj = sum(tmp.*p,'all') + sum(c.*y,'all');

constraints = [A'*p + c*diag(lambda) == v, A*y <= tmp*diag(lambda), lambda(:)>=0, p(:)>=0];    

optimize(constraints, obj, options);

obj = value(obj)/nsamp;