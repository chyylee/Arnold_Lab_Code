%% Parallel_RESAMPLE_EN_multinomial.m
% GOAL:To perform Resampling of data and fix imbalance in class distribution.

% INPUT: Load Input from make_xyblocks.m.
% This will give:
% patients
% xblock
% yblock (as a #patients x 2 logical)

% REQUIRES: re_glmnet_multi_en.m

% OUTPUT: Excel file with sorted selected features by frequency selected as
% well as a workspace with the results.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Lee, Arnold Lab, University of Michigan, Biomedical Engineering
% April 22nd, 2020
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% LOAD DATA AND DECIDE NUM RUNS AND ALPHA
% evrimovepath('bottom') % must move PLS Toolbox to bottom of path
clear;
filename = 'For_2016_2020-06-21_COVID&Family_posvsFamily_negvsHealthyNo_Mod__log_PCA-PLSDA.mat';
load(filename)
num_runs = 400; % How many times to repeate resampling
alpha_val = 1; % alpha value

% Determines number of classes in input
num_classes = size(yblock,2);

% Selects smallest class size as resampling size:
[min_size, idx] = min(sum(yblock));
max_size = max(sum(yblock));
ssize = min_size;

% k-fold is adjusted by sample size: MAY NEED TO ADJUST 
kfold_val = floor(2*ssize/5); % Number of partitions
inner_kf = floor(2*ssize/5); % Inner K-fold validation

tot_resam_runs = ceil(num_runs/kfold_val);

seed_val = 1:tot_resam_runs; 



%% Split observations for the classes and resample at the same size (size 
% of smaller class)
clear C_dat
for x_id = 1:num_classes
    temp = xblock(yblock(:,x_id),:);
    if size(temp,1) < max_size
        dif = max_size - size(temp,1);
        temp = [temp; NaN(dif,size(temp,2))];
    end
    C_dat(x_id,:,:) = temp;
end

% Generate the resampled datasets
all_C_dat = zeros(length(seed_val),num_classes,ssize,size(xblock,2));

for i = 1:length(seed_val)
    rng(seed_val(i));
    for cl_id = 1:num_classes
        temp = squeeze(C_dat(cl_id,:,:));
        if sum(sum(isnan(temp))) ~= 0
            tnan = isnan(squeeze(C_dat(cl_id,:,:)));
            temp(tnan(:,1),:) = [];
        end
        [sm_dat, sm_dat_idx] = datasample(temp,ssize,'replace',false);
        all_C_dat(i,cl_id,:,:) = sm_dat;
    end

end

% Prepare empty matrices to for parallel computing of the lasso results
all_B_1SE = zeros(length(seed_val),kfold_val,size(xblock,2));
all_B_min = zeros(length(seed_val),kfold_val,size(xblock,2));

fc_count = [0 0 0 0]; % Count warning labels

%% Uncomment to force a predictor into a model
% idx_npen = 1; % Makes sure the feature at this index is included in model
% pf = ones(size(xblock,2),1);
% pf(idx_npen) = 1;
%% Run all combinations of resampled data on glmnet()/LASSO-EN
for yset = 1:num_classes
    combine_y(yset*ssize - (ssize -1):yset*ssize) = yset*ones(ssize,1);
    temp_x(:,yset*ssize - (ssize -1):yset*ssize,:) = squeeze(all_C_dat(:,yset,:,:));
end

myPool = parpool('local',6); 

parfor set_id = 1:length(seed_val)
    combine_x = squeeze(temp_x(set_id,:,:));
    rng(set_id) 
    
    c = cvpartition(combine_y,'Kfold',kfold_val);
    for part_id = 1:kfold_val
        temp_idx = training(c,part_id);
        xblock = combine_x(temp_idx,:);
        yblock = combine_y(temp_idx)';
        
        % Call glmnet to run LASSO/EN
        [fc,B_1SE,B_min,~,~] = re_glmnet_multi_en(xblock,yblock,xnames,alpha_val,inner_kf,[]);
        all_B_1SE(set_id,part_id,:) = B_1SE;
        all_B_min(set_id,part_id,:) = B_min;
        
        fc_count = fc_count + fc; % Keeps track if run had non-empty selection
        
        [set_id part_id] % Counter for runs
    end
end

delete(myPool)

%% Reshapes LASSO/EN results and and evaluate

c = 1;
re_B_1SE = zeros(size(all_B_1SE,1)*size(all_B_1SE,2),size(all_B_1SE,3));
re_B_min = zeros(size(all_B_1SE,1)*size(all_B_1SE,2),size(all_B_1SE,3));
for i = 1:size(all_B_1SE,1)
    for j = 1:size(all_B_1SE,2)
        re_B_1SE(c,:) = all_B_1SE(i,j,:);
        re_B_min(c,:) = all_B_min(i,j,:);
        c = c + 1;
    end
end

% Trim unwanted data
if size(re_B_1SE,1) ~= num_runs
    re_B_1SE(num_runs+1:end,:) = [];
    re_B_min(num_runs+1:end,:) = [];
end

% Calculates selection frequency from total runs
% fr_sel_1SE = sum(re_B_1SE ~= 0)./(size(all_B_1SE,1)*size(all_B_1SE,2));
% fr_sel_min = sum(re_B_min ~= 0)./(size(all_B_1SE,1)*size(all_B_1SE,2));


% Calculates selection frequency from only successful runs
fr_sel_1SE = sum(re_B_1SE ~= 0)./(fc_count(1) + fc_count(3));
fr_sel_min = sum(re_B_min ~= 0)./(fc_count(1) + fc_count(2));



%% Compile Data and get ready to save
T = array2table([fr_sel_1SE', fr_sel_min'],'RowNames',cellstr(xnames), ...
    'VariableNames',{'EN_1SE','EN_min'});

t1 = T(T.EN_1SE >= 0.7,1);
t2 = T(T.EN_min >= 0.7,1);

T_sorted_1SE = sortrows(T,1,'descend');
T_sorted_min = sortrows(T,2,'descend');

tot_mdls = size(all_B_1SE,1)*size(all_B_1SE,2);

%% Visualize Results
num_sel = 15;
idxF = T.EN_1SE > 0.7; % Greater than 70% Selection

new_var = exp(re_B_min(:,idxF));
freq_var = T.EN_min(idxF);

figure()
subplot(1,2,1)
boxplot(new_var, 'orientation', 'horizontal');
set(gca,'ytick',1:length(idxF),'yticklabels',xnames(idxF));
hold on
% vline(1,'k')
% vl = xline(0);
line([1 1], [0 length(idxF)]);
vl.Color = 'k';
xlabel('Odds Ratio')
title(strcat('Coefficients for top ', ' ',num2str(num_sel), ' ', ' Features in', num2str(tot_mdls),' Mdls'))

subplot(1,2,2)
b = barh(freq_var);
b.FaceColor = [0.75 0.75 0.75];
set(gca,'ytick',1:length(idxF),'yticklabels',xnames(idxF));
hold on
%vline(0.7,'k')
line([1 1], [0 0.7]);
xlabel('Frequency of Selection')
title(strcat('Frequency for top',  ' ', num2str(num_sel), ' ', ' Features in ', num2str(tot_mdls),' Mdls'))

