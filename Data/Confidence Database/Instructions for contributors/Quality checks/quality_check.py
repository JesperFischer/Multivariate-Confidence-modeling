# -*- coding: utf-8 -*-
"""
--------------------------------------------------------------------------
 Quality check for contributions to the Confidence Database

 Please, run this file as a basic quality check for your contribution(s).
 Note that this file will not catch all potential errors and you should
 carefully inspect your contribution before sending it.

 Instructions:
    - Place this file and all files you're sending in a new folder.
    - Run this file. If there are any errors or warnings, update your
    files until you have addressed all issues.
    - Send your contribution(s) and note that you ran this quality check.

 Original MATLAB Code written by Doby Rahnev.
 Translated to Python 3 by Alan Lee. Last update: Nov 14, 2019.
 Versions of python/libraries:
 - python: 3.5.6
 - pandas: 0.23.4
--------------------------------------------------------------------------
"""

import glob
import pandas as pd

''' Get all data files from database folder '''

data_files = sorted([f for f in glob.glob('data_*.csv')])
readme_files = sorted([f for f in glob.glob('readme_*.txt')])


''' Load spreadsheet '''
T = pd.read_excel('Spreadsheet_template.xlsx')
n_datasets = len(T)
names_in_spreadsheet = sorted(list(T.Name_in_database))


''' Determine if equal number of files are present in the spreadsheet, data and readme '''

num_names = len(names_in_spreadsheet)
num_data = len(data_files)
num_readme = len(readme_files)
numList = [num_names, num_data, num_readme]
print('# entres in spreadsheet: %d' % num_names)
print('# data files: %d' % num_data)
print('# readme files: %d' % num_readme)
if not(all([n==numList[0] for n in numList])):
    raise('ERROR. Number of files doesn''t match!\n')
else:
    print('OK: Number of files matches.\n')


''' Determine if all names match '''

for dataset_num, data_file in enumerate(data_files):
    spreadsheetname = names_in_spreadsheet[dataset_num]
    datafilename = data_file[5:-4]
    readmename = readme_files[dataset_num][7:-4]
    match_spreadsheet_datafile = spreadsheetname == datafilename
    match_spreadsheet_readme = spreadsheetname == readmename
    if not(match_spreadsheet_datafile) or not(match_spreadsheet_readme):
        print('Name in spreadsheet: %s' % spreadsheetname)
        print('Name in data files: %s' % datafilename)
        print('Name in readme files: %s' % readmename)
        raise('ERROR. Name doesn''t match!\n')
print('OK: Names are consistent between the files and the spreadsheet.\n');


''' Determine if data are loading well and if fields are named correctly '''

print('\n-----Checking if the data are loading well and if data columns have correct names.------');
pd.set_option('display.max_columns', 500) # force pandas to show 500 columns
for dataset_num, datafile in enumerate(data_files):

    # load the dataset
    dataframe = pd.read_csv(datafile)

    # Display how the data is being read in
    print('\nDataset name: %s' % names_in_spreadsheet[dataset_num])
    print('This is how the data from this dataset are being read in:')
    print(dataframe.head(8)) # show first 8 rows
    print('Please check that columns that should be numeric are indeed numeric.')
    print('If not, this most likely indicates a problem with the formatting of the data.')

    # Check if different fields are present and give warnings if not
    warningMessage = '\nWARNING! No field %s exists. Please use the exact word %s unless you have to use something else.\n'
    if not('Subj_idx' in dataframe.keys()):
        raise('ERROR. No field "Subj_idx" exists. A field "Subj_idx" MUST be present in the dataset.')
    if not('Stimulus' in dataframe.keys()):
        print(warningMessage % ('"Stimulus"','"Stimulus"'));
    if not('Response' in dataframe.keys()):
        print(warningMessage % ('"Response"','"Response"'));
    if not('Confidence' in dataframe.keys()):
        print(warningMessage % ('"Confidence"','"Confidence"'));
    if not('RT_dec' in dataframe.keys()) and not('RT_decConf' in dataframe.keys()):
        print(warningMessage % ('"RT_dec" or "RT_decConf"','"RT_dec" or "RT_decConf"'));
    if not('RT_conf' in dataframe.keys()) and 'RT_dec' in dataframe.keys():
        print(warningMessage % ('"RT_conf"','"RT_conf"'));

print('')

''' Determine if number of subjects are reported correctly in spreadsheet '''

print('\n-----Checking if reported number of subjects is correct.------');
for dataset_num, datafile in enumerate(data_files):

    # load the dataset
    dataframe = pd.read_csv(datafile)

    # get number of subjects from different sources
    numsubj_spreadsheet = T.Num_subjects[dataset_num]
    numsubj_datafile = len(set(dataframe.Subj_idx))

    # Display basic info
    print('Dataset name: %s' % names_in_spreadsheet[dataset_num])
    print('Number of subjects reported in spreadsheet: %d' % numsubj_spreadsheet)
    print('Number of actual subjects in data: %d\n' % numsubj_datafile)

    # Check for inconsistency in number of subjects
    if numsubj_spreadsheet != numsubj_datafile:
        raise('ERROR. Number of subjects doesn''t match between spreadsheet and actual data.')


''' Determine if number of trials per subject are reported correctly in spreadsheet '''

print('\n-----Checking if reported number of trials per subject is correct.------');
for dataset_num, datafile in enumerate(data_files):

    # load the dataset
    dataframe = pd.read_csv(datafile)

    # Determine the number of trials per subject
    subject_names = sorted(list(set(dataframe.Subj_idx)))
    trials_per_subj = [
            len(dataframe.Response[dataframe.Subj_idx==subjname]) for subjname in subject_names]

    # Display basic info
    print('Dataset name: %s' % names_in_spreadsheet[dataset_num])
    print('Min total trials per subject reported in spreadsheet: %d' % T.Min_trials_per_subject[dataset_num])
    print('Min total trials per subject in the actual data: %d\n' % min(trials_per_subj))
    print('Max total trials per subject reported in spreadsheet: %d' % T.Max_trials_per_subject[dataset_num])
    print('Max total trials per subject in the actual data: %d\n' % max(trials_per_subj))

    # Check for inconsistency in number of subjects
    if T.Min_trials_per_subject[dataset_num] != min(trials_per_subj):
        raise('ERROR. The min total numbers per subject doesn''t match between spreadsheet and actual data.')
    if T.Max_trials_per_subject[dataset_num] != max(trials_per_subj):
        raise('ERROR. The max total numbers per subject doesn''t match between spreadsheet and actual data.')
