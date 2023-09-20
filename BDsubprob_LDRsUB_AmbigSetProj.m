function [obj, rtime, mu, Omega] = BDsubprob_LDRsUB_AmbigSetProj(A, b, c, gamma1, gamma2)

k = size(A,1);
Omega = sdpvar(k,k);
M = sdpvar(k+1,k+1);
mu = sdpvar(k,1);

obj = trace(A*Omega) + b'*mu + c;
constraints = [norm(mu) <= sqrt(gamma1), gamma2*eye(k)-Omega >= 0];
constraints = [constraints, M == [1 mu'; mu Omega]];
constraints = [constraints, M >= 0];
options = sdpsettings('verbose',0,'solver','mosek','cachesolvers',1);

optimize(constraints, -obj, options); % maximization so '-obj'
rtime = ans.solvertime;  
obj = value(obj);
mu = value(mu);
Omega = value(Omega);

