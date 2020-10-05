%% CREATE HEATMAP FOR FEATURE SELECTED VARIABLES
% May-16-2020

%%
clear;
% LOAD FEATURE SELECTED WORKSPACE

load('21-Jun-2020smaple_PLSDA_Model.mat')
% LOAD ORIGINAL DATA FROM MAKE XY
orig_data = load('2020-06-21_COVID&Family_posvsFamily_negvsHealthyNo_Mod__log_PCA-PLSDA.mat');
classes = orig_data.classes;
patients = orig_data.patients;
yblock = orig_data.yblock;

% LOAD DESIRED COLORMAPS
load('colormaps.mat')

% FORMAT VARIABLES FOR CLUSTERGRAM
% alter if you only want a subset of patients or features
 index_patients = 1:length(orig_data.patients);%yblock(:, 1);
 index_features = 1:length(orig_data.xnames);
 
class_vector = patients(index_patients); %just to get the right size
 for k = 1:length(classes)
     class_vector(yblock(index_patients,k)==1) = classes{k};
 end
 
data = xblock;

%% CLUSTERGRAM: COLOR COLUMNS BY AGE COHORT

% Define Colors for Age Cohort
num_class = size(yblock,2); % Number of classes
Y = zeros(size(yblock,1),1);
for cat_num = 1:num_class
    Y = Y + yblock(:,cat_num).*cat_num;
end

% cmap = [1 1 1; 0.75 0.75 1; 0.2 0.2 0.2];
% color_define = cell(size(Y));
% color_define(Y >= 30) = {cmap(1,:)}; % ELDERLY
% color_define(Y < 30) = {cmap(2,:)}; % Adults
% % color_define(Y == 3) = {cmap(3,:)}; % Children

% Cont. Y variable (not categorical)
cmap = gray(length(Y));
color_define = cell(size(Y));
for i = 1:length(Y)
    color_define(i) = {cmap(i,:)};
end

s = struct('Labels',cellstr(num2str(Y)),'Colors',color_define);
pedit = strrep(patients,'.','_');
%s = struct('Labels',cellstr(patients),'Colors',color_define);

normdata = zscore(data)';


% Clustergram
c1 = clustergram(normdata, 'RowLabels', cellstr(xnames), ...
    'ColumnLabels',Y,...
    'Colormap',myredblue,'Cluster','all');


% set(gcf, 'Renderer', 'Painters');

c1.LabelsWithMarkers = true;
c1.ColumnLabelsColor = s;

%% CLUSTERGRAM: COLOR COLUMN lABLES BY GENDER
% 
% % LOAD ANALAGOUS WS TO GET GENDER INFORMATION
% gend_dat = load('16-Apr-2020_No-tet-flu_FvsM_ADULTS&ELDERLY&KIDS_log_PCA-PLSDA.mat');
% y = gend_dat.yblock;
% 
% % Define Colors for Female/Male
% cmap = [1 0 1; 0 0 1];
% Ygender = y(:,1);
% color_define = cell(size(Ygender));
% color_define(Ygender == 1) = {cmap(1,:)}; % Female
% color_define(Ygender == 0) = {cmap(2,:)}; % Male
% sg = struct('Labels',cellstr(string(Ygender)),'Colors',color_define);
% 
% % Clustergram
% c1 = clustergram(normdata, 'RowLabels', cellstr(xnames), ...
%     'ColumnLabels',cellstr(string(Ygender)),...
%     'Colormap',myredblue,'Cluster','all');
% 
% c1.LabelsWithMarkers = true;
% c1.ColumnLabelsColor = sg;

