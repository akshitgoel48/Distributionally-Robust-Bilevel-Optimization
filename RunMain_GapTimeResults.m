clc;
clear all;
warning('off','all')

% ------------------------ Specify Parameters ----------------------------

inst = {}; i=0;
for d = [15, 20, 25, 30, 35]    % no. of locations considered    
    for cost = [750, 1000]
        if d <= 30
            i = i+1;
            inst{i} = [d, 5, 2, cost];
            i = i+1;
            inst{i} = [d, 5, 3, cost];   
        end
        if d >= 25
            i = i+1;
            inst{i} = [d, 10, 5, cost];        
        end 
    end        
end
inst = inst([1:8,11,14]);

%
base_folder = strcat(pwd,'/computational time_and_gap');
if not(isfolder(base_folder))
    mkdir(base_folder);
end
base_folder = strcat(base_folder, '/');
nsamp = 10;

%
C_Vneq0 = 0;
GenResults_Gap_Time(C_Vneq0, base_folder, nsamp, inst)

%
C_Vneq0 = 1;
GenResults_Gap_Time(C_Vneq0, base_folder, nsamp, inst)
