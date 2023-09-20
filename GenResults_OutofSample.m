function GenResults_OutofSample(d, s, t, r, base_folder, C_Vneq0)

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
end

%%
[~, k] = size(W);
T2 = 10; % # randomizations out-of-sample data
N2 = 5000; % out-of-sample data size

xi_uniform = zeros(k, N2,T2); % of size (k,N2)
xi_uniform2 = zeros(k, N2,T2); % of size (k,N2)
xi_normal = zeros(k, N2,T2); % of size (k,N2)

[Lb2, Ub2] = cal_support_unif(0.92*mu_true, 0.80*var_true);

cl = parcluster('local');
pool = cl.parpool(min(cl.NumWorkers, T2));
parfor t2 = 1:T2
    rng(t2);
    xi_out_uniform = Lb + (Ub - Lb)*rand([k,N2]); 
    xi_out_uniform2 = Lb2 + (Ub2 - Lb2)*rand([k,N2]); 
    xi_out_normal = zeros(k,N2);
    for num = 1:N2
    xi_out_normal(:, num) = mu_true + sqrt(var_true)*trandn((Lb-mu_true)/sqrt(var_true)*ones(k,1),... 
                            (Ub-mu_true)/sqrt(var_true)*ones(k,1));
    end    
    
    xi_uniform(:,:,t2) = xi_out_uniform;        % of size (k,N2) 
    xi_uniform2(:,:,t2) = xi_out_uniform2;
    xi_normal(:,:,t2) = xi_out_normal;          % of size (k,N2) 
end
delete(pool)

solution_set = {};
solution_set{1} = zeros(d,1);
poss_indices = find(s-t==1);    index = 1;

for loop1=1:size(poss_indices,1)
    temp_indices = nchoosek(poss_indices,loop1);
    for loop2 = 1:size(temp_indices,1)
        temp = zeros(d,1);
        temp(temp_indices(loop2,:)) = 1;
        index = index + 1;
        solution_set{index} = temp;
    end
end

%
second_stage_cost_uniform = {};
second_stage_cost_uniform2 = {};
second_stage_cost_normal = {}; 

cl = parcluster('local');
pool = cl.parpool(min(cl.NumWorkers, T2));
for index = 1:size(solution_set,2)       
    obj_normal_temp = zeros(T2,1);
    obj_uniform_temp = zeros(T2,1);
    obj_uniform2_temp = zeros(T2,1);
    
    xopt_LDR = solution_set{index};
    
    parfor t2 = 1:T2
        xi_out_uniform = xi_uniform(:,:,t2); % of size (k,N2) 
        xi_out_uniform2 = xi_uniform2(:,:,t2); % of size (k,N2) 
        xi_out_normal = xi_normal(:,:,t2); % of size (k,N2)  
        
        obj_uniform_temp(t2) =  second_stage_2(xopt_LDR, xi_out_uniform, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0);
        obj_uniform2_temp(t2) =  second_stage_2(xopt_LDR, xi_out_uniform2, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0);
        obj_normal_temp(t2) = second_stage_2(xopt_LDR, xi_out_normal, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0);   
    end
    second_stage_cost_uniform{index} = obj_uniform_temp;  %#ok<*AGROW,*SAGROW>
    second_stage_cost_uniform2{index} = obj_uniform2_temp;
    second_stage_cost_normal{index} = obj_normal_temp;
end
delete(pool)

T1 = 100; % # randomizations in-sample data
N1mat = [10, 100, 1000]; % matrix of in-sample data sample sizes

Mean_TSBP_1 = zeros(T1,size(N1mat,2));
Mean_TSBP_2 = zeros(T1,size(N1mat,2));
Mean_TSBP_3 = zeros(T1,size(N1mat,2));

BOX_TSBP_1 = zeros(T1*T2,size(N1mat,2));
BOX_TSBP_2 = zeros(T1*T2,size(N1mat,2));
BOX_TSBP_3 = zeros(T1*T2,size(N1mat,2));

Mean_DRBP1_1 = zeros(T1,size(N1mat,2));
Mean_DRBP1_2 = zeros(T1,size(N1mat,2));
Mean_DRBP1_3 = zeros(T1,size(N1mat,2));

BOX_DRBP1_1 = zeros(T1*T2,size(N1mat,2));
BOX_DRBP1_2 = zeros(T1*T2,size(N1mat,2));
BOX_DRBP1_3 = zeros(T1*T2,size(N1mat,2));

Mean_DRBP2_1 = zeros(T1,size(N1mat,2));
Mean_DRBP2_2 = zeros(T1,size(N1mat,2));
Mean_DRBP2_3 = zeros(T1,size(N1mat,2));

BOX_DRBP2_1 = zeros(T1*T2,size(N1mat,2));
BOX_DRBP2_2 = zeros(T1*T2,size(N1mat,2));
BOX_DRBP2_3 = zeros(T1*T2,size(N1mat,2));

if C_Vneq0 == 1
    base_folder = strcat(base_folder, 'C_Vneq0_out_samp_obj');
else
    base_folder = strcat(base_folder, 'C_Veq0_out_samp_obj');
end

if not(isfolder(base_folder))
    mkdir(base_folder);
end
path = strcat(base_folder,'/');

%
cl = parcluster('local');
pool = cl.parpool(min(cl.NumWorkers, T1));
for index = 1:size(N1mat,2)

N1 = N1mat(index);

disp(strcat('Out-of-sample Tests Running for in-sample size = ', num2str(N1)))

[act_obj_uniform, act_obj_uniform2, act_obj_normal] = out_of_sample_obj(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                                                        Lb, Ub, NaN, NaN, 'TS-BP', N1, T1, T2, solution_set, second_stage_cost_uniform, second_stage_cost_uniform2, second_stage_cost_normal);

Mean_TSBP_1(:,index) = mean(act_obj_uniform,2);
Mean_TSBP_2(:,index) = mean(act_obj_uniform2,2);
Mean_TSBP_3(:,index) = mean(act_obj_normal,2);

BOX_TSBP_1(:,index) = reshape(act_obj_uniform',[T1*T2,1]);
BOX_TSBP_2(:,index) = reshape(act_obj_uniform2',[T1*T2,1]);
BOX_TSBP_3(:,index) = reshape(act_obj_normal',[T1*T2,1]);

[act_obj_uniform, act_obj_uniform2, act_obj_normal] = out_of_sample_obj(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                                                        Lb, Ub, 0, 1, 'DRBP', N1, T1, T2, solution_set, second_stage_cost_uniform, second_stage_cost_uniform2, second_stage_cost_normal);

Mean_DRBP1_1(:,index) = mean(act_obj_uniform,2);
Mean_DRBP1_2(:,index) = mean(act_obj_uniform2,2);
Mean_DRBP1_3(:,index) = mean(act_obj_normal,2);

BOX_DRBP1_1(:,index) = reshape(act_obj_uniform',[T1*T2,1]);
BOX_DRBP1_2(:,index) = reshape(act_obj_uniform2',[T1*T2,1]);
BOX_DRBP1_3(:,index) = reshape(act_obj_normal',[T1*T2,1]);

[act_obj_uniform, act_obj_uniform2, act_obj_normal] = out_of_sample_obj(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                                                        Lb, Ub, 0.5, 1, 'DRBP', N1, T1, T2, solution_set, second_stage_cost_uniform, second_stage_cost_uniform2, second_stage_cost_normal);

Mean_DRBP2_1(:,index) = mean(act_obj_uniform,2);
Mean_DRBP2_2(:,index) = mean(act_obj_uniform2,2);
Mean_DRBP2_3(:,index) = mean(act_obj_normal,2);

BOX_DRBP2_1(:,index) = reshape(act_obj_uniform',[T1*T2,1]);
BOX_DRBP2_2(:,index) = reshape(act_obj_uniform2',[T1*T2,1]);
BOX_DRBP2_3(:,index) = reshape(act_obj_normal',[T1*T2,1]);

disp(' ')
end
delete(pool)

%
writematrix(N1mat, strcat(path,'N1mat.csv'));

writematrix(Mean_TSBP_1, strcat(path,'Mean_TSBP_uniform.csv'));
writematrix(Mean_TSBP_2, strcat(path,'Mean_TSBP_uniform2.csv'));
writematrix(Mean_TSBP_3, strcat(path,'Mean_TSBP_normal.csv'));

writematrix(BOX_TSBP_1, strcat(path,'BOX_TSBP_uniform.csv'));
writematrix(BOX_TSBP_2, strcat(path,'BOX_TSBP_uniform2.csv'));
writematrix(BOX_TSBP_3, strcat(path,'BOX_TSBP_normal.csv'));

writematrix(Mean_DRBP1_1, strcat(path,'Mean_DRBP1_uniform.csv'));
writematrix(Mean_DRBP1_2, strcat(path,'Mean_DRBP1_uniform2.csv'));
writematrix(Mean_DRBP1_3, strcat(path,'Mean_DRBP1_normal.csv'));

writematrix(BOX_DRBP1_1, strcat(path,'BOX_DRBP1_uniform.csv'));
writematrix(BOX_DRBP1_2, strcat(path,'BOX_DRBP1_uniform2.csv'));
writematrix(BOX_DRBP1_3, strcat(path,'BOX_DRBP1_normal.csv'));

writematrix(Mean_DRBP2_1, strcat(path,'Mean_DRBP2_uniform.csv'));
writematrix(Mean_DRBP2_2, strcat(path,'Mean_DRBP2_uniform2.csv'));
writematrix(Mean_DRBP2_3, strcat(path,'Mean_DRBP2_normal.csv'));

writematrix(BOX_DRBP2_1, strcat(path,'BOX_DRBP2_uniform.csv'));
writematrix(BOX_DRBP2_2, strcat(path,'BOX_DRBP2_uniform2.csv'));
writematrix(BOX_DRBP2_3, strcat(path,'BOX_DRBP2_normal.csv'));