%% Convert scripts for MATLAB2016
% MATLAB 2016 cannot handle string information so need to convert to cell
% and then load
clear;
%% Load desired workspace

ws_name = '2020-06-21_COVID&Family_posvsFamily_negvsHealthyNo_Mod__log_PCA-PLSDA.mat';
load(ws_name)

% A = whos;
% for i = 1:length(A)
%     A(1).class
% end

xnames = cellstr(xnames);
patients = cellstr(patients);
ttl = cellstr(ttl);
filettl = cellstr(filettl);

save(strcat('For_2016_',ws_name))