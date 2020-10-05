function [fc_count, B_1SE, B_min,xnames_1SE, xnames_min] = re_glmnet_multi_en(xblock,yblock,xnames,ALPHA,innerkf,pf)
%% Pre-process Y-block
% Convert the 2 column logical to one column

num_class = size(yblock,2); % Number of classes
Y = zeros(size(yblock,1),1);
for cat_num = 1:num_class
    Y = Y + yblock(:,cat_num).*cat_num;
end

% Define class names
% class_names = classes;

%% Call multinomial glmnet

options=glmnetSet;

% Set Alpha
options.alpha=ALPHA; % Elastic Net (0 = Ridge, 1 = LASSO)

% Define which variable to modify penalty;
% idx_npen = 1; % Makes sure the feature at this index is included in model
% pf = ones(size(xblock,2),1);
% pf(idx_npen) = 1;
pf = ones(size(xnames));
options.penalty_factor=pf;

% Mean Center and Variance Scale
Xnorm = zscore(xblock);

% call glmnet for a 2-way classifier (binomial), 'auc' as error metric
fit=cvglmnet(Xnorm,Y,'multinomial',options,'class',innerkf); 

% Diagnostic Plot
cvglmnetPlot(fit); 

% Best model is where the CVE is the lowest. Looks through fit for the min CVE and records the index of the row (index) where CVE is lowest
[bestcvm,idx_min]=min(fit.cvm); 
numbas_features=fit.glmnet_fit.df(idx_min); %goes to row of lowest MSE and pulls out features

% model_min = fit.glmnet_fit.beta(:,idx_min);

% 1SE from the lowest CVE (tries to reduce number of freatures while staying within 1SE of min CVE)
[idx_1se] = find(fit.lambda_1se == fit.glmnet_fit.lambda); 
num_1se = fit.glmnet_fit.df(idx_1se);

% model_1SE =fit.glmnet_fit.df(:,idx_1se);


%%
beta_lambda = fit.glmnet_fit.beta;
best_model_min = zeros(size(beta_lambda{:,1},1),size(beta_lambda,2));
for i = 1:size(best_model_min,2)
    temp = beta_lambda{i};
    best_model_min(:,i) = temp(:,idx_min);
end

beta_lambda = fit.glmnet_fit.beta;
best_model_1se = zeros(size(beta_lambda{:,1},1),size(beta_lambda,2));
for i = 1:size(best_model_min,2)
    temp = beta_lambda{i};
    best_model_1se(:,i) = temp(:,idx_1se);
end
% B_min = fit.glmnet_fit.beta(:,idx_min);
% B_1SE = fit.glmnet_fit.beta(:,idx_1se);

[r_min c_min] = find(abs(best_model_min)>0); %Finds all the nonzero values in the best model

% Delineates which class each coefficient is non-zero
tot_feat = zeros(1,num_class);
for i = 1:num_class
    tot_feat(i) = sum(c_min == i);
end
tot_feat_per_class = tot_feat;

% Define new xblock
xblock_min = xblock(:,unique(r_min));
xnames_min = xnames(unique(r_min));

[r_1se c_1se] = find(abs(best_model_1se)>0); %Finds all the nonzero values in the best model

% Delineates which class each coefficient is non-zero
tot_feat = zeros(1,num_class);
for i = 1:num_class
    tot_feat(i) = sum(c_1se == i);
end
tot_feat_per_class = tot_feat;

% Define new xblock
xblock_1SE = xblock(:,unique(r_1se));
xnames_1SE = xnames(unique(r_1se));
%%
B_1SE = sum(abs(best_model_1se),2);
B_min = sum(abs(best_model_min),2);
% % Chart with Beta Values

Beta_vals = array2table([B_min B_1SE],'rownames',cellstr(xnames),...
    'variablenames',{'B_min','B_1SE'});

loc_min_feat = B_min ~= 0;
loc_1SE_feat = B_1SE ~= 0;

%%

if size(xblock_min,2) == 0 && size(xblock_1SE,2) == 0
    WFLAG = 'Warning: Both min MSE and 1SE signatures contain no predictors';
    disp(WFLAG)
    fc_count = [0 0 0 1];
elseif size(xblock_1SE,2) == 0 && size(xblock_min,2) ~= 0
    WFLAG = 'Warning: 1SE signature contains no predictors';
    disp(WFLAG)
    fc_count = [0 1 0 0];
elseif size(xblock_1SE,2) ~= 0 && size(xblock_min,2) == 0
    WFLAG = 'Warning: min MSE signature contains no predictors';
    disp(WFLAG)
    fc_count = [0 0 1 0];
else
    WFLAG = 'Successful selection';
    disp(WFLAG);
    fc_count = [1 0 0 0];
    disp(strcat(num2str(sum(B_min ~= 0)),' feature signature'))
end