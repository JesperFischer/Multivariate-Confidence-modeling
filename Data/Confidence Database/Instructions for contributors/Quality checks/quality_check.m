%--------------------------------------------------------------------------
% Quality check for contributions to the Confidence Database
%
% Please, run this file as a basic quality check for your contribution(s).
% Note that this file will not catch all potential errors and you should
% carefully inspect your contribution before sending it.
%
% Instructions:
%    - Place this file and all files you're sending in a new folder.
%    - Run this file. If there are any errors or warnings, update your
%    files until you have addressed all issues.
%    - Send your contribution(s) and note that you ran this quality check.
%
% Code written by Doby Rahnev. Last update: Nov 14, 2019.
%--------------------------------------------------------------------------

clear
clc

% Get all data files from database folder
data_files = dir(fullfile(pwd, 'data_*.csv'));
readme_files = dir(fullfile(pwd, 'readme_*.txt'));

% Load spreadsheet
T = readtable('Spreadsheet_template.xlsx');
n_datasets = length(T.Name_in_database);
names_in_spreadsheet = sort(T.Name_in_database);


%% Determine if equal number of files are present in the spredsheet, data and readme
fprintf(['# entries in spreadsheet: ' num2str(length(names_in_spreadsheet)) ...
    '\n# data files: ' num2str(length(data_files)) ...
    '\n# readme files: ' num2str(length(readme_files)) '\n']);
if length(data_files) ~= length(names_in_spreadsheet) || length(data_files) ~= length(readme_files)
    fprintf('ERROR. Number of files doesn''t match!\n');
    return
end


%% Determine if all names match
for dataset_num=1:length(data_files)
    if ~strcmp(names_in_spreadsheet{dataset_num}, data_files(dataset_num).name(6:end-4)) || ...
            ~strcmp(names_in_spreadsheet{dataset_num}, readme_files(dataset_num).name(8:end-4))
        fprintf(['\nERROR. Name doesn''t match!' ...
            '\nName in spreadsheet: ' names_in_spreadsheet{dataset_num} ...
            '\nName in data files: ' data_files(dataset_num).name(6:end-4) ...
            '\nName in readme files: ' readme_files(dataset_num).name(8:end-4) '\n']);
        return
    end
end
fprintf('\nNames are consistent between the files and the spreadsheet.\n');


%% Determine if data are loading well and if fields are named correctly
fprintf('\n-----Checking if the data are loading well and if data columns have correct names.------\n');
for dataset_num=1:length(data_files)
    
    % Load the dataset
    table = readtable(fullfile(pwd, data_files(dataset_num).name));
    
    % Display how the data is being read in
    fprintf(['\nDataset name: ' names_in_spreadsheet{dataset_num} ...
        '\nThis is how the data from this dataset are being read in:\n']);
    head(table)
    fprintf(['\nPlease check that columns that should be numeric are indeed numeric. If not, '...
        'this most likely indicates a problem with the formatting of the data.\n']);
    
    % Check if different fields are present and give warnings if not
    if ~ismember('Subj_idx', table.Properties.VariableNames)
        fprintf('ERROR. No field "Subj_idx" exists. A field "Subj_idx" MUST be present in the dataset.\n');
        return
    end
    if ~ismember('Stimulus', table.Properties.VariableNames)
        fprintf('\n\nWARNING! No field "Stimulus" exists. Please use the exact word "Stimulus" unless you have to use something else.\n\n');
    end
    if ~ismember('Response', table.Properties.VariableNames)
        fprintf('\n\nWARNING! No field "Response" exists. Please use the exact word "Response" unless you have to use something else.\n\n');
    end
    if ~ismember('Confidence', table.Properties.VariableNames)
        fprintf('\n\nWARNING! No field "Confidence" exists. Please use the exact word "Confidence" unless you have to use something else.\n\n');
    end
    if ~ismember('RT_dec', table.Properties.VariableNames) && ~ismember('RT_decConf', table.Properties.VariableNames)
        fprintf('\n\nWARNING! No field "RT_dec" or "RT_decConf" exists. Please use the exact word "RT_dec" or "RT_decConf" unless you have to use something else.\n\n');
    end
    if ~ismember('RT_conf', table.Properties.VariableNames) && ismember('RT_dec', table.Properties.VariableNames)
        fprintf('\n\nWARNING! No field "RT_conf" exists. Please use the exact word "RT_conf" unless you have to use something else.\n\n');
    end
end


%% Determine if number of subjects are reported correctly in spreadsheet
fprintf('\n-----Checking if reported number of subjects is correct.------\n');
for dataset_num=1:length(data_files)
    
    % Load the dataset
    table = readtable(fullfile(pwd, data_files(dataset_num).name));
    
    % Display basic info
    fprintf(['Dataset name: ' names_in_spreadsheet{dataset_num} ...
        '\nNumber of subjects reported in spreadsheet: ' num2str(T.Num_subjects(dataset_num)) ...
        '\nNumber of actual subjects in data: ' num2str(length(unique(table.Subj_idx))) '\n\n']);
    
    % Check for inconsistencies in number of subjects
    if T.Num_subjects(dataset_num) ~= length(unique(table.Subj_idx))
        fprintf('ERROR. Number of subjects doesn''t match between spreadsheet and actual data.\n');
        return
    end
end


%% Determine if number of trials per subject are reported correctly in spreadsheet
fprintf('-----Checking if reported number of trials per subject is correct.------\n');
for dataset_num=1:length(data_files)
    
    % Load the dataset
    clear trials_per_subj
    table = readtable(fullfile(pwd, data_files(dataset_num).name));
    
    % Determine the number of trials per subject
    subject_names = unique(table.Subj_idx);
    for sub=1:length(unique(table.Subj_idx))
        if isnumeric(subject_names(sub))
            trials_per_subj(sub) = sum(table.Subj_idx==subject_names(sub));
        else
            trials_per_subj(sub) = sum(strcmp(table.Subj_idx,subject_names(sub)));
        end
    end
    
    % Display basic info
    fprintf(['Dataset name: ' names_in_spreadsheet{dataset_num} ...
        '\nMin total trials per subject reported in spreadsheet: ' num2str(T.Min_trials_per_subject(dataset_num)) ...
        '\nMin total trials per subject in the actual data: ' num2str(min(trials_per_subj)) ...
        '\nMax total trials per subject reported in spreadsheet: ' num2str(T.Max_trials_per_subject(dataset_num)) ...
        '\nMax total trials per subject in the actual data: ' num2str(max(trials_per_subj)) '\n\n']);
    
    % Check for inconsistencies in number of subjects
    if T.Min_trials_per_subject(dataset_num) ~= min(trials_per_subj)
        fprintf('ERROR. The min total numbers per subject doesn''t match between spreadsheet and actual data.\n');
        return
    end
    if T.Max_trials_per_subject(dataset_num) ~= max(trials_per_subj)
        fprintf('ERROR. The max total numbers per subject doesn''t match between spreadsheet and actual data.\n');
        return
    end
end