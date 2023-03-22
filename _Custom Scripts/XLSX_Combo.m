%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                XLSX Combo                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <r.akhundov@griffith.edu.au>

%Custom script for the EMG_Classifier: Read a participant-folder containing 
%classification results, combine all, save results in new combo_xlsx file
%in the previously chose folder.

%%%Requirements: 
%1)Classification results for all participants completed
%  (c3d files should have participant number/identifier in their 
%   name for clarity)

clc; close all; clearvars;
try
%% Get Data
fsp = filesep;
folder = uigetdir(path,'Select Acquisition Folder');
fileList = dir([folder, fsp,'**',fsp,'*.xlsx']);

%Check if selected folder is correct
if isempty(fileList) == 1
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%% ERROR: No .xlsx files in selected folder %%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    return
end

%Change the current folder to the folder of this m-file (if user does "add to path" instead of "Change folder")
if(~isdeployed) 
  cd(fileparts(which(mfilename)));
end

tic;
disp('%% Started combining classification results %%');disp('%')
%Turn unnecessary warnings off 
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
warning('off', 'MATLAB:DELETE:FileNotFound')

%XLSX directories
[~,acquisition] = fileparts(folder); %Acquisition name
old_XLSX = [pwd, fsp, 'Classifications_Template.xlsx'];
new_XLSX = [folder, fsp, 'Combined_Classifications_', acquisition, '.xlsx'];
delete(new_XLSX); %Delete old xlsx version to overwrite it
copyfile (old_XLSX,new_XLSX);

%% Read Data
XLSX_tables = cell(length(fileList),1);
for s=1:length(fileList)
var_tables = readtable(fullfile(fileList(s).folder, fileList(s).name),'Sheet','Sheet1');
XLSX_tables{s,1} = table2cell(var_tables(1:end-4,:));
end

%% Check for equal size & combine
%Check
if length(unique(cellfun('size',XLSX_tables,2))) ~= 1
   %Helper variable to find inconsistencies 
   incons_helper = cell(length(fileList)+1,2);
   incons_helper(1,1) = cellstr('Participant XLSX'); incons_helper(1,2) = cellstr('Number of EMGs');
   incons_helper(2:end,1) = {fileList.name}'; 
   incons_helper(2:end,2) = num2cell(cellfun('size',XLSX_tables,2));
   
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
   disp('%%%% ERROR: Inconsistent number of EMGs between participants. %%%%')
   disp('%%%%        Check the incons_helper variable.                 %%%%')
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
   return
end

%Combine
XLSX_tables = vertcat(XLSX_tables{:});

%% Make cell for new xlsx and save
header_XLSX_tables = var_tables.Properties.VariableNames; %Muscle names

%Adjust muscle names if readtable changed it
if header_XLSX_tables{2}(1) == 'x'
    header_XLSX_tables(2:end) = cellfun(@(x) x(2:end), header_XLSX_tables(2:end), 'un', 0);
end

%Combine and save
combo_XLSX_tables = vertcat(header_XLSX_tables, XLSX_tables);
xlswrite(new_XLSX,combo_XLSX_tables,'Sheet1','A1');

catch err
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%% ERROR while creating new.xlsx %%%%')
        disp(err.message)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
t1 = toc;   
disp(['%% XLSX_Combo: Finished everything successfully in ' num2str(t1) ' seconds %%'])
