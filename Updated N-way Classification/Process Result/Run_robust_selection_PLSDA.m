%% robust_selection.m
% GOAL: Pick the number of features to include in the model based on the 
% frequency at which the features are observed by resampling.

% INPUT: Load results from resampling (ie Parallel_Resample.m). Make sure
% to feed in the original xblock before any feature selection, original
% xnames and original yblock.

% REQUIRES: PLSDA_TB.m, PCA_TB.m ScoresandLoadingsPlot.m,
% colorschemes.mat

% OUTPUT: The best number of features to select based on resampling.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Lee, Arnold Lab, University of Michigan, Biomedical Engineering
% April 26th, 2020
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Load the Resampling Output
clear;
%evrimovepath('top')
fn = 'Resample_Output.mat';
load(fn)

%% Rerun plsda with different combinations of highest selected features

% Only test for 50% of total features
p_trys = 0.5;
num_trys = floor(size(xblock,2)*p_trys);
% Get list of best features based on 1SE frequency
[bfeat, bidx] = maxk(fr_sel_1SE,num_trys);

class_train = NaN(num_trys,1);
class_test = NaN(num_trys,1);

for feat_id = 1:num_trys
    sel_idx = bidx(1:feat_id);
    X = xblock(:,sel_idx);
    Y = yblock;
    sel_nms = xnames(sel_idx);
    [tr,te] = PLSDA_TB(X,Y,xnames(sel_idx),fn,classes,yname,ttl,0);
    close(gcf)
    close(gcf)
    if feat_id == 1
        class_train(feat_id) = tr(1);
        class_test(feat_id) = te(1);
    else
        class_train(feat_id) = tr(2);
        class_test(feat_id) = te(2);
    end

    disp(['Iteration Number: ', num2str(feat_id)])
end

%% Find best result (step above zero)
[minval, minidx] = mink(class_test,num_trys);

dif_val = class_test - class_test(minidx(1));

thresh = 0.1;
best_sel = find(dif_val < thresh, 1, 'first');

disp(['Best Number of Features = ', num2str(best_sel)])
%% Find best result (step above zero)
[minval, minidx] = mink(class_test,num_trys);

dif_val = class_test - class_test(minidx(1));

thresh = 0.1;
best_sel = find(dif_val < thresh, 1, 'first');
if best_sel==1
    best_sel = 1+find(dif_val(2:end) < thresh, 1, 'first');
end
while dif_val(best_sel)==dif_val(best_sel+1)
best_sel=best_sel+1;
end


 
%% Plot Result
feat_num = 1:num_trys;

p1 = plot(feat_num,class_train,'-+');
hold on
p2 = plot(feat_num,class_test,'-+');
p3 = plot(best_sel,class_train(best_sel),'s','MarkerEdgeColor','k', ...
    'MarkerSize',10);
p4 = plot(best_sel,class_test(best_sel),'s','MarkerEdgeColor','k', ...
    'MarkerSize',10);
legend([p1, p2,p3],'CalErr','CVErr','Selected')
xlabel('Number of Features')
ylabel('Error')
set(gca,'fontsize',14)

% SAVE SELECTED MODEL
sel_idx = bidx(1:best_sel);
X = xblock(:,sel_idx);
Y = yblock;
nw_xnames = xnames(sel_idx);
%%
wsfilettitle = strcat(datestr(today()),strtok(extractAfter(filettl,11),'.'),'_ROB_SEL_MDL.mat');
save(wsfilettitle)

%% Call Melissa's PLSDA Code

[tr,te] = PLSDA_TB(X,Y,xnames(sel_idx),fn,classes,yname,ttl,1);

%% Call PCA Code
fn = '2020-06-16-all_4_categories_resampling_result.mat'
PCA_TB(xblock,yblock,xnames,fn,classes,yname,ttl)