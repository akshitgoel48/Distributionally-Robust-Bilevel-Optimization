clc;
clear all;
warning('off','all')

% ------------------------ Specify Parameters ----------------------------

d = 8;   % no. of locations considered    
   
loc = 5;    % no. of potential locations for retail stores

locA = 2;    % no. of retail stores operated by A  

s = zeros([d,1]);
s([1,2,3,4,6]) = 1;   % choose potential locations for retail stores

t = zeros([d,1]);
t(6) = 1;    % choose which retail stores operated by A

r = sum(s)-sum(t);   % max no. of stores B can open

%
base_folder = strcat(pwd,'/Boxplots');
if not(isfolder(base_folder))
    mkdir(base_folder);
end
base_folder = strcat(base_folder, '/');

C_Vneq0 = 0;
GenResults_OutofSample(d, s, t, r, base_folder, C_Vneq0)

%
C_Vneq0 = 1;
GenResults_OutofSample(d, s, t, r, base_folder, C_Vneq0)
