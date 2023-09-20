function [Lb1, Ub1] = cal_support_unif(mu1, var1)

for i=1:size(mu1,2)
    syms u l;
    eqns = [u - l == sqrt(12*var1(i)), u + l == 2*mu1(i)];
    S=solve(eqns, [u,l]);
    Lb1(i) = double(S.l);
    Ub1(i) = double(S.u);
end