function GenResults_Gap_Time(C_Vneq0, base_folder, nsamp, inst)

date_time = datetime('now','Format','d-MMM-y-HH_mm_ss');

if C_Vneq0 == 1
    base_folder = strcat(base_folder, 'C_Vneq0_');
else
    base_folder = strcat(base_folder, 'C_Veq0_');
end

gamma_mat = {};
gamma_mat{1} = [0, 1]; gamma_mat{2} = [1, 1]; 

alpha = 1; % === alpha === 

folder_name = strcat(base_folder, string(date_time),'/');

for index = 1:size(gamma_mat,2)  % === gamma1 === , === gamma2 === 
    
gamma1 = gamma_mat{index}(1);  gamma2 = gamma_mat{index}(2);  

folder = strcat(folder_name,'GAMMA1_',num2str(gamma1),'_GAMMA2_',num2str(gamma2),'/');
if not(isfolder(folder))
    mkdir(folder);
end

st_gap2 = struct('SNo', {}, 'd', {}, 's', {}, 't', {}, 'r', {}, 'w_i', {}, ...
                'SDP25', {}, 'SDP50', {}, 'SDP75', {}, 'NULL1',{}, ...
                'PROJ25', {}, 'PROJ50', {}, 'PROJ75', {}, 'NULL2',{}, ...
                'IA_COP25', {}, 'IA_COP50', {}, 'IA_COP75', {});
            
st_gap3 = struct('SNo', {}, 'd', {}, 's', {}, 't', {}, 'r', {}, 'w_i', {}, ...              
                'SDP_MPtime', {}, 'SDP_SPtime', {}, 'SDPtime', {}, 'SDP_Iters', {}, 'SDP_AvgGap', {}, 'NULL1',{}, ...
                'PROJ_MPtime', {}, 'PROJ_SPtime', {}, 'PROJtime', {}, 'PROJ_Iters', {}, 'PROJ_AvgGap', {}, 'NULL2',{}, ...
                'IA_COP_MPtime', {},'IA_COP_SPtime', {}, 'IA_COPtime', {}, 'IA_COP_Iters', {}, 'IA_COP_AvgGap', {}, 'NULL3',{}, ...
                'DiscreteLB_MPtime', {},'DiscreteLB_SPtime', {}, 'DiscreteLB_time', {}, ...
                'DiscreteLB_Iters', {});            
    
for num = 1:size(inst,2)  % === num === 

vec = inst{num};

d = vec(1); 

loc = vec(2);    % no. of potential locations for retail stores

locA = vec(3);    % no. of retail stores operated by A  

cost = vec(4);  % coefficient of x

if C_Vneq0 == 1
    cost = 2*cost;
end

st_gap = struct('SNo', {}, 'd', {}, 's', {}, 't', {}, 'r', {}, 'w_i', {}, ...
                'Discrete_LB', {}, 'Discrete_LB_gap', {}, 'DiscreteLB_Iters', {}, ...
                'DiscreteLB_MPtime',{},'DiscreteLB_SPtime',{},'DiscreteLB_rtime', {}, 'NULL1',{}, ...
                'SDP_xsol', {}, 'SDP', {}, 'SDP_gap', {}, 'SDP_Iters', {}, ...
                'SDP_MPtime',{},'SDP_SPtime',{}, 'SDP_rtime', {}, 'Gap_SDP', {}, 'NULL2',{}, ...
                'PROJ_xsol', {}, 'PROJ', {}, 'PROJ_gap', {}, 'PROJ_Iters', {}, ...
                'PROJ_MPtime',{},'PROJ_SPtime',{}, 'PROJ_rtime', {}, 'Gap_PROJ', {}, 'NULL3',{}, ...
                'IA_COP_xsol', {}, 'IA_COP', {}, 'IA_COP_gap', {}, 'IA_COP_Iters', {}, ...
                'IA_COP_MPtime',{},'IA_COP_SPtime',{}, 'IA_COP_rtime', {}, 'Gap_IA_COP', {});

Lb = 50.0;
Ub = 150.0;

% true distribution: uniform 
mu_true = 0.5*(Lb+Ub);
var_true = (Ub-Lb)^2/12;

set_size = 10;      % number of random instances for each setting

st_gap(set_size).d = d;
st_gap(set_size).s = loc;
st_gap(set_size).t = locA;
st_gap(set_size).r = loc - locA;
st_gap(set_size).w_i = cost;

st_gap2(num).d = d;
st_gap2(num).s = loc;
st_gap2(num).t = locA;
st_gap2(num).r = loc - locA;
st_gap2(num).w_i = cost;

st_gap3(num).d = d;
st_gap3(num).s = loc;
st_gap3(num).t = locA;
st_gap3(num).r = loc - locA;
st_gap3(num).w_i = cost;

cl = parcluster('local');
pool = cl.parpool(min(cl.NumWorkers, set_size));
parfor sets = 1:set_size  % === sets ===     
        
    rng(sets);
    ind = randperm(d,loc);   % 'randomly' choose potential locations for retail stores

    s = zeros([d,1]);
    s(ind) = 1;
    
    t = zeros([d,1]);
    t(ind(randperm(loc,locA))) = 1;    % 'randomly' choose which retail stores operated by A
        
    r = sum(s)-sum(t);   % max no. of stores B can open
    
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

    % generate cost c for the follower's problem
    Coordinates = rand([2,d]);
    c0 = zeros([I*J,1]);
    for i = 1:I
        for j = 1:J
            c0((i-1)*J+j) = norm(Coordinates(:,listI(i)) - Coordinates(:,listJ(j)));
        end
    end
            
    % cost v for the leader
    v0 = -5*ones([I*J,1]);
    
    bigM = 1e6;    

    capA = Ub*(d-sum(s))/sum(t); 
    capB = alpha*capA; 

    b20 = zeros([I,1]);
    for i = 1:I
       if t(listI(i)) == 1
           b20(i) = capA;
       end
    end

    b2i = cell([d,1]);
    position = 1;
    for i = 1:d
        b2i{i} = zeros([I,1]);
        if s(i) == 1
            b2i{i}(position) = capB; 
            position = position + 1;
        end
    end

    b0 = [b10; b20; b30]; 

    % support W\xi >= h
    W = [eye(J); -eye(J)];
    h = [Lb*ones([J,1]); -Ub*ones([J,1])];        
             
    C = zeros(I*J, J); V = zeros(I*J, J); w = cost*ones(d,1);
    if C_Vneq0 == 1
        C = 0.1*(0.9 + 0.02*rand(I*J,J))/(Ub*J);
        V = -20*rand(I*J, J)/(Ub*J);
    end
    
    for i = 1:I
        if t(listI(i)) == 1 
            for j = 1:J
                v0((i-1)*J+j) = 0;       
                V((i-1)*J+j,:) = 0;
            end
        end
    end

    %----------------------------------------------------------------------
    disp(['================*Currently running set ',num2str(sets),'*================'])
    
    st_gap(sets).SNo = sets;  
    st_gap(sets).Gap_SDP = 0;
    st_gap(sets).Gap_PROJ = 0;
    st_gap(sets).Gap_IA_COP = 0;
    
    num_Iters = 1;
    
    for iters = 1:num_Iters    % === numIters ===
    rng(iters);
    xi = Lb+(Ub-Lb)*rand(J,10);
    Mean1 = mean(xi,2);     % calculate mean
    Cov1 = cov(xi');     % calculate covariance   
       
    %--------------------------------------------------------------------------------------------------------------------------------------------
    
    [obj, x3, rtime, gap_BD, master_time, subprob_time, numIters] = BDmaster_LDRsUB_SDPbyS_lemma(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                                                                 gamma1, gamma2, Mean1, Cov1, bigM, 1e-4, Inf, 'gurobi', 0);   
    st_gap(sets).SDP_xsol = find(round(x3)==1);             st_gap(sets).SDP = obj;              st_gap(sets).SDP_rtime = rtime;
    st_gap(sets).SDP_gap = gap_BD;
    st_gap(sets).SDP_MPtime = master_time;      st_gap(sets).SDP_SPtime = subprob_time; 
    st_gap(sets).SDP_Iters = numIters;
    disp(['*Ended running DRBP_SDP for set ',num2str(sets),'* *Time Elapsed*: ',num2str(rtime),', *Obj*: ', num2str(obj)])
    
    %--------------------------------------------------------------------------------------------------------------------------------------------
    [obj, x3, rtime, gap_BD, master_time, subprob_time, numIters] = BDmaster_LDRsUB_AmbigSetProj(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                                                                                 gamma1, gamma2, Mean1, Cov1, bigM, 0);
    st_gap(sets).PROJ_xsol = find(round(x3)==1);             st_gap(sets).PROJ = obj;              st_gap(sets).PROJ_rtime = rtime;
    st_gap(sets).PROJ_gap = gap_BD;
    st_gap(sets).PROJ_MPtime = master_time;      st_gap(sets).PROJ_SPtime = subprob_time; 
    st_gap(sets).PROJ_Iters = numIters;
    disp(['*Ended running DRBP_PROJ for set ',num2str(sets),'* *Time Elapsed*: ',num2str(rtime),', *Obj*: ', num2str(obj)])                                     
    
    %--------------------------------------------------------------------------------------------------------------------------------------------
    [obj, x3, rtime, gap_IA_COP_BD, master_time, subprob_time, numIters] = BDmaster_LDRsUB_IA_COP(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                                                                  gamma1, gamma2, Mean1, Cov1, bigM);
    st_gap(sets).IA_COP_xsol = find(round(x3)==1);             st_gap(sets).IA_COP = obj;              st_gap(sets).IA_COP_rtime = rtime;
    st_gap(sets).IA_COP_gap = gap_IA_COP_BD;
    st_gap(sets).IA_COP_MPtime = master_time;       st_gap(sets).IA_COP_SPtime = subprob_time; 
    st_gap(sets).IA_COP_Iters = numIters;
    disp(['*Ended running DRBP_IA_COP for set ',num2str(sets),'* *Time Elapsed*: ',num2str(rtime),', *Obj*: ', num2str(obj)])
          
    %--------------------------------------------------------------------------------------------------------------------------------------------
    
    numIters_LB = 1; 
    st_gap(sets).Discrete_LB = -Inf;
    st_gap(sets).DiscreteLB_rtime = 0;
    st_gap(sets).Discrete_LB_gap = 0;
    
    for itersLB = 1:numIters_LB    % === numIters_LB ===  
    rng(itersLB);    
    xi_supp_base = Lb + (Ub - Lb)*rand([J,nsamp]);
    xi_supp = [xi, xi_supp_base];    % Discrete Support set for lower bound computation.
    
    [obj_val, ~, rtime, gap_BD, master_time, subprob_time, numIters] = BDmaster_DiscrtLB_SDP(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W,... 
                                                                                             gamma1, gamma2, xi, xi_supp, bigM, 0, 1);
    
    disp(['*Ended computing Disceret_LB_BD',num2str(itersLB),' for set ',num2str(sets),'* *Time Elapsed*: ',num2str(rtime),', *Obj*: ', num2str(obj_val)])
        
    if obj_val > st_gap(sets).Discrete_LB
        st_gap(sets).Discrete_LB = obj_val;
    end 
    
    st_gap(sets).DiscreteLB_rtime = st_gap(sets).DiscreteLB_rtime + rtime;
    st_gap(sets).Discrete_LB_gap = st_gap(sets).Discrete_LB_gap + gap_BD;
    st_gap(sets).DiscreteLB_MPtime = master_time; st_gap(sets).DiscreteLB_SPtime = subprob_time; 
    st_gap(sets).DiscreteLB_Iters = numIters;
    end  % === numIters_LB ===  
    
    st_gap(sets).DiscreteLB_rtime = st_gap(sets).DiscreteLB_rtime/numIters_LB;
    st_gap(sets).Discrete_LB_gap = st_gap(sets).Discrete_LB_gap/numIters_LB;    
           
    if round(st_gap(sets).Discrete_LB) == 0 && round(st_gap(sets).SDP) == 0
        st_gap(sets).Gap_SDP = st_gap(sets).Gap_SDP + 0;
    else
        st_gap(sets).Gap_SDP = st_gap(sets).Gap_SDP + (-st_gap(sets).Discrete_LB + st_gap(sets).SDP)/abs(st_gap(sets).Discrete_LB)*100;
    end   
    
    if round(st_gap(sets).Discrete_LB) == 0 && round(st_gap(sets).PROJ) == 0
        st_gap(sets).Gap_PROJ = st_gap(sets).Gap_PROJ + 0;
    else
        st_gap(sets).Gap_PROJ = st_gap(sets).Gap_PROJ + (-st_gap(sets).Discrete_LB + st_gap(sets).PROJ)/abs(st_gap(sets).Discrete_LB)*100;
    end  
    
    if round(st_gap(sets).Discrete_LB) == 0 && round(st_gap(sets).IA_COP) == 0
        st_gap(sets).Gap_IA_COP = st_gap(sets).Gap_IA_COP + 0;
    else
        st_gap(sets).Gap_IA_COP = st_gap(sets).Gap_IA_COP + (-st_gap(sets).Discrete_LB + st_gap(sets).IA_COP)/abs(st_gap(sets).Discrete_LB)*100;
    end   
    end  % === numIters ===
    
    st_gap(sets).Gap_SDP = st_gap(sets).Gap_SDP/num_Iters;
    st_gap(sets).Gap_PROJ = st_gap(sets).Gap_PROJ/num_Iters;
    st_gap(sets).Gap_IA_COP = st_gap(sets).Gap_IA_COP/num_Iters;
    
end  % === sets === 
delete(pool)
Tb_gap = struct2table(st_gap);
writetable(Tb_gap, strcat(folder,'Gap_Results',num2str(num),'.csv'));

st_gap2(num).SDP25 = quantile(Tb_gap.Gap_SDP, 0.25); 
st_gap2(num).SDP50 = quantile(Tb_gap.Gap_SDP, 0.50); 
st_gap2(num).SDP75 = quantile(Tb_gap.Gap_SDP, 0.75);
st_gap3(num).SDP_MPtime = mean(Tb_gap.SDP_MPtime);   
st_gap3(num).SDP_SPtime = mean(Tb_gap.SDP_SPtime); 
st_gap3(num).SDPtime = mean(Tb_gap.SDP_rtime);
st_gap3(num).SDP_Iters = mean(Tb_gap.SDP_Iters); 
st_gap3(num).SDP_AvgGap = mean(Tb_gap.Gap_SDP);

st_gap2(num).PROJ25 = quantile(Tb_gap.Gap_PROJ, 0.25); 
st_gap2(num).PROJ50 = quantile(Tb_gap.Gap_PROJ, 0.50); 
st_gap2(num).PROJ75 = quantile(Tb_gap.Gap_PROJ, 0.75);
st_gap3(num).PROJ_MPtime = mean(Tb_gap.PROJ_MPtime);   
st_gap3(num).PROJ_SPtime = mean(Tb_gap.PROJ_SPtime); 
st_gap3(num).PROJtime = mean(Tb_gap.PROJ_rtime);
st_gap3(num).PROJ_Iters = mean(Tb_gap.PROJ_Iters); 
st_gap3(num).PROJ_AvgGap = mean(Tb_gap.Gap_PROJ); 

st_gap2(num).IA_COP25 = quantile(Tb_gap.Gap_IA_COP, 0.25); 
st_gap2(num).IA_COP50 = quantile(Tb_gap.Gap_IA_COP, 0.50); 
st_gap2(num).IA_COP75 = quantile(Tb_gap.Gap_IA_COP, 0.75);
st_gap3(num).IA_COP_MPtime = mean(Tb_gap.IA_COP_MPtime);   
st_gap3(num).IA_COP_SPtime = mean(Tb_gap.IA_COP_SPtime); 
st_gap3(num).IA_COPtime = mean(Tb_gap.IA_COP_rtime);
st_gap3(num).IA_COP_Iters = mean(Tb_gap.IA_COP_Iters); 
st_gap3(num).IA_COP_AvgGap = mean(Tb_gap.Gap_IA_COP);

st_gap3(num).DiscreteLB_MPtime = mean(Tb_gap.DiscreteLB_MPtime);   
st_gap3(num).DiscreteLB_SPtime = mean(Tb_gap.DiscreteLB_SPtime); 
st_gap3(num).DiscreteLB_time = mean(Tb_gap.DiscreteLB_rtime); 
st_gap3(num).DiscreteLB_Iters = mean(Tb_gap.DiscreteLB_Iters); 

if num >= 2
    Tb_gap2 = struct2table(st_gap2);
    writetable(Tb_gap2, strcat(folder,'Gap_QuantileResults.csv'));
    Tb_gap3 = struct2table(st_gap3);
    writetable(Tb_gap3, strcat(folder,'TimeResults.csv'));
end
end  % === num === 

end  % === gamma1 === , === gamma2 === 
