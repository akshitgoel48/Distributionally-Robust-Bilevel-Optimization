function L_xi = lift_space(xi, Mean1)
J = size(Mean1,1);
L_xi = zeros(2*J,size(xi,2));

for j = 1:size(xi,2)
    for i=1:size(xi,1)
        L_xi(2*i-1,j) = min(xi(i,j),Mean1(i));
        L_xi(2*i,j) = max(xi(i,j)-Mean1(i),0);
    end
end