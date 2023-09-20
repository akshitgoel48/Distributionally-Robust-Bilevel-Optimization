function GenResults_InSample(d, s, t, r, base_folder, C_Vneq0)

cost = 950; 
w = cost* ones([d,1]); % coefficient of x

J = d - sum(s);  % number of demand locations
I = sum(s);  % number of potential stores
listI = zeros([I,1]); % list of potential stores
listJ = zeros([J,1]); % list of demand locations

count1 = 1;
count2 = 1;
for i = 1:d
    if s(i) == 1
        listI(count1) = i;
        count1 = count1 + 1;
    else
        listJ(count2) = i;
        count2 = count2 + 1;
    end
end

% matrices and vectors for constraint sum_i y_{ij} >= \xi_j
A1 = zeros([J,I*J]);
for i = 1:J
    for count = 0:I-1
        A1(i,i+count*J) = -1.0;
    end
end

B1i = zeros(J);
B10 = -eye(J);
b1i = zeros([J,1]);
b10 = zeros([J,1]);

% matrices and vectors for constraint sum_j y_{ij} <= b_i(x_i + t_i)
A2 = zeros([I,I*J]);
for i = 1:I
    for count = 1:J
        A2(i,(i-1)*J+count) = 1.0;
    end
end
B2i = zeros([I,J]);
B20 = zeros([I,J]);

% matrices and vectors for y >= 0
A3 = -eye(I*J);
B3i = zeros([I*J,J]);
B30 = zeros([I*J,J]);
b3i = zeros([I*J,1]);
b30 = zeros([I*J,1]);

A = [A1; A2; A3];
B0 = [B10; B20; B30];

c0 = readmatrix('c0.csv');

% cost v for the leader
v0 = -5*ones([I*J,1]);

bigM = 1e6;
  
Lb = 30.0;
Ub = 240.0;

C = readmatrix('C.csv');
V = readmatrix('V.csv');
for i = 1:I
    if t(listI(i)) == 1 
        for j = 1:J
            v0((i-1)*J+j) = 0;
            V((i-1)*J+j,:) = 0;
        end
    end
end

gamma1_mat = [0 0.2 0.5]; 
gamma2_mat = [1 3]; 

bi = Ub*(d-sum(s))/sum(t); 

b20 = zeros([I,1]);
for i = 1:I
   if t(listI(i)) == 1
       b20(i) = bi;
   end
end

b2i = cell([d,1]);
position = 1;
for i = 1:d
    b2i{i} = zeros([I,1]);
    if s(i) == 1
        b2i{i}(position) = bi/2; % bi;
        position = position + 1;
    end
end

b0 = [b10; b20; b30]; 

% support W\xi >= h
W = [eye(J); -eye(J)];
h = [Lb*ones([J,1]); -Ub*ones([J,1])];

% true distribution: uniform 
mu_true = 0.5*(Lb+Ub);
var_true = (Ub-Lb)^2/12;

Mean1 = mu_true*ones(J,1);  % calculate mean
Cov1 = var_true*eye(J); % calculate covariance

if C_Vneq0 == 0
    cost = 305;
    C = zeros(size(C)); V = zeros(size(V)); w = cost*ones(d,1);
    gamma1_mat = [0 1 1.5]; 
    gamma2_mat = [1 2];   
end

% 
disp('================*Currently running in-sample codes*================')

xi = readmatrix('ksi_in-sample.csv'); 
Mean2 = mean(xi,2);     % calculate mean
Cov2 = cov(xi');     % calculate covariance

% --------------------------SIMPLE TS-BP----------------------------------
                        
[obj, x5, ~, ~, rtime] = Yanikoglu_Kuhn_PessTSBP(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h, Lb, Ub,... 
                                                 Mean2, Cov2, bigM, 0, 0);   

disp(['*Ended running Yanikoglu_Kuhn_PessTSBP, *Time Elapsed*: ',num2str(rtime),', *Obj*: ', num2str(obj)])
disp('Leader opens facility at locations:')
disp(find(round(x5)==1)');
disp(' ');

for gamma2 = gamma2_mat

for gamma1 = gamma1_mat

    disp(['================*GAMMA1*: ',num2str(gamma1),', *GAMMA2*: ',num2str(gamma2),'*================'])  
        
    [obj, x3, rtime, ~] = BDmaster_LDRsUB_SDPbyS_lemma(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                       gamma1, gamma2, Mean2, Cov2, bigM, 1e-4, Inf, 'gurobi', 0);     
    disp(['*Ended running SDPbyS_lemma, *Time Elapsed*: ',num2str(rtime),', *Obj*: ', num2str(obj)]) % ,...
    disp('Leader opens facility at locations:')
    disp(find(round(x3)==1)');        
    disp(' ');
    
    [obj, x3, rtime, ~, ~, ~, Iters] = BDmaster_LDRsUB_AmbigSetProj(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                                                    gamma1, gamma2, Mean2, Cov2, bigM, 0);
    
    disp(['*Ended running DRBP_Proj, *Time Elapsed*: ',num2str(rtime),', *Iters*: ',num2str(Iters),', *Obj*: ', num2str(obj)]) % ,...
    disp('Leader opens facility at locations:')
    disp(find(round(x3)==1)');        
    disp(' ');
end
end

if C_Vneq0 == 0
    folder = strcat(base_folder, 'sensitivity_plots_CVeq0');
else
    folder = strcat(base_folder, 'sensitivity_plots_CVneq0');
end
if not(isfolder(folder))
    mkdir(folder);
end

if C_Vneq0 == 0
    % C = V = 0
    gamma1_mat = 1.1:0.1:1.5;
    gamma2_mat = 1:0.25:2.0;
else
    % C, V neq 0
    gamma1_mat = 0.2:0.05:0.5;
    gamma2_mat = 1:0.25:10;
end

%% SENSITIVITY ANALYSIS w.r.t. gamma1 & gamma2
insamp_objval = zeros(size(gamma1_mat,2),size(gamma2_mat,2));
    
xi_in = readmatrix('ksi_in-sample.csv');  
Mean1 = mean(xi_in,2);     % calculate mean
Cov1 = cov(xi_in');     % calculate covariance

for i1 = 1:size(gamma1_mat,2)
for i2 = 1:size(gamma2_mat,2)
    gamma1 = gamma1_mat(i1);
    gamma2 = gamma2_mat(i2);
    [insamp_objval(i1,i2), ~, ~, ~] = BDmaster_LDRsUB_SDPbyS_lemma(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                                       gamma1, gamma2, Mean1, Cov1, bigM, 1e-4, Inf, 'gurobi', 0); 
end
end

% ---------------------------------------------------------------------------
if C_Vneq0 == 0
    % C = V = 0
    writematrix(gamma1_mat, strcat(folder,'/C_Veq0_gamma1_mat.csv'))
    writematrix(gamma2_mat, strcat(folder,'/C_Veq0_gamma2_mat.csv'))
    writematrix(insamp_objval, strcat(folder,'/C_Veq0_insamp_obj.csv'))
    gamma2_mat_lb = 1:0.25:2.0;
    gamma1 = 1.3;
else
    % C, V neq 0
    writematrix(gamma1_mat, strcat(folder,'/C_Vneq0_gamma1_mat.csv'))
    writematrix(gamma2_mat, strcat(folder,'/C_Vneq0_gamma2_mat.csv'))
    writematrix(insamp_objval, strcat(folder,'/C_Vneq0_insamp_obj.csv'))
    gamma2_mat_lb = [1, 3, 4, 5, 6, 7];
    gamma1 = 0.2;
end

%% SENSITIVITY to support's lower bound
column = 0;
Lb1 = [0,  10, 20, 30, 40, 50];
Ub1 = Ub*ones(size(Lb1));

insamp_objval = zeros(size(Lb1,2),size(gamma2_mat_lb,2));

xi_in = readmatrix('ksi_in-sample.csv');  
Mean1 = mean(xi_in,2);     % calculate mean
Cov1 = cov(xi_in');     % calculate covariance

cl = parcluster('local');
pool = cl.parpool(min(cl.NumWorkers, size(Lb1,2)));
for gamma2 = gamma2_mat_lb
column = column+1;
parfor i = 1:size(Lb1,2)
	h1 = [Lb1(i)*ones(J,1); -Ub1(i)*ones(J,1)];
    [obj, ~, ~, ~] = BDmaster_LDRsUB_SDPbyS_lemma(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h1,...
                                                  gamma1, gamma2, Mean1, Cov1, bigM, 1e-4, Inf, 'gurobi', 0);         
    insamp_objval(i, column) = obj;  
end
end
delete(pool)

if C_Vneq0 == 0
    writematrix(gamma2_mat_lb, strcat(folder,'/C_Veq0_gamma2_mat_lb.csv'))
    writematrix(Lb1, strcat(folder,'/C_Veq0_Lb_mat.csv'))
    writematrix(insamp_objval, strcat(folder,'/C_Veq0_gamma1_',num2str(gamma1),'_insamp_obj.csv'))
else
    writematrix(gamma2_mat_lb, strcat(folder,'/C_Vneq0_gamma2_mat_lb.csv'))
    writematrix(Lb1, strcat(folder,'/C_Vneq0_Lb_mat.csv'))
    writematrix(insamp_objval, strcat(folder,'/C_Vneq0_gamma1_',num2str(gamma1),'_insamp_obj.csv'))
end
