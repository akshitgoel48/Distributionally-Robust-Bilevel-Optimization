function [act_obj_uniform, act_obj_uniform2, act_obj_normal] = out_of_sample_obj(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,... 
                                            Lb, Ub, gamma1, gamma2, model, N1, T1, T2, solution_set, second_stage_cost_uniform, ...
                                            second_stage_cost_uniform2, second_stage_cost_normal)
                     
bigM = 1e6;

[~, k] = size(W);
num_solution = size(solution_set,2);

act_obj_normal = zeros(T1,T2);
act_obj_uniform = zeros(T1,T2);
act_obj_uniform2 = zeros(T1,T2);

parfor t1 = 1:T1
    rng(t1);
    xi_in = Lb + (Ub - Lb)*rand([k,N1]);  % of size (k,N1)
    Mean1 = mean(xi_in,2);     % calculate mean
    Cov1 = cov(xi_in');     % calculate covariance
    
    xopt_LDR = zeros(d,1);
    if strcmp(model,'TS-BP')
        disp(strcat('Running Yanikoglu_Kuhn_PessTSBP Rnadom instance No. ', num2str(t1)));
        [~, xopt_LDR, ~, ~, ~] = Yanikoglu_Kuhn_PessTSBP(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h, Lb, xi_in, Mean1, Cov1, bigM, 0, 0);        

    elseif strcmp(model,'DRBP')
        disp(strcat('Running SDPbyS_lemma Random instance No. ', num2str(t1)));
        [~, xopt_LDR, ~, ~] = BDmaster_LDRsUB_SDPbyS_lemma(d, w, s, t, r, A, b1i, b2i, b3i, B0, b0, v0, C, V, c0, W, h,...
                                                           gamma1, gamma2, Mean1, Cov1, bigM, 1e-4, Inf, 'gurobi', 0); 
    end
    
    xopt_LDR = round(xopt_LDR);
    
    for index = 1:num_solution
        if isequal(xopt_LDR, solution_set{index})    %#ok<*PFBNS>
            break;            
        end
    end        

    act_obj_uniform(t1,:) =  w'*xopt_LDR + second_stage_cost_uniform{index}';
    act_obj_uniform2(t1,:) =  w'*xopt_LDR + second_stage_cost_uniform2{index}';
    act_obj_normal(t1,:) = w'*xopt_LDR + second_stage_cost_normal{index}';
end    


