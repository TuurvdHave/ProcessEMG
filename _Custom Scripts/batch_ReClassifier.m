%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            batch_ReClassifier                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <r.akhundov@griffith.edu.au>

%Custom script for the EMG_Classifier: Read new location of image files and 
%write into adjusted XLSX for all participants.

%%%Requirements: 
%1)Classification results for participants completed 
%  (EMG_Classifier successfully executed)

clc; close all; clearvars;
try
%% Get Data from Participant Folders
fsp = filesep;
dirParticipant = uigetdir(path,'Select Participant Folder Containing .c3d Files');

XLSXdir = dir([dirParticipant, fsp,'**',fsp,'*.xlsx']);

%Check if selected folder is correct
if isempty(XLSXdir) == 1
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%% ERROR: No .xlsx file in selected folder %%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    return
end

tic;
for s=1:length(XLSXdir)
dirClassifier = XLSXdir(s).folder(1:end-5);
disp(['%% Started adjusting classification results #', num2str(s), ' %%']);disp('%')
%Turn unnecessary warnings off 
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
warning('off', 'MATLAB:DELETE:FileNotFound')
warning('off', 'MATLAB:DELETE:Permission')

%Get directories into workspace
dirImages = [dirClassifier, fsp, 'Images'];
dirXLSX = [dirClassifier, fsp, 'XLSX'];

%Get participant
[~,participant] = fileparts(dirClassifier);
participant = strsplit(participant,'_');

%XLSX directories
old_XLSX = [dirXLSX, fsp,'Classifications_', participant{2}, '.xlsx'];
new_XLSX = [dirXLSX, fsp,'Adjusted_Classifications_', participant{2}, '.xlsx'];
copyfile (old_XLSX,new_XLSX);

%Get images
imds = imageDatastore(dirImages, 'IncludeSubfolders', true, 'LabelSource', 'foldernames');

%% New Data Matrix Creation for XLSX
new_class_matrix = cell(length(imds.Files),3);

%Get image names
for i=1:length(imds.Files)
    [~,img_name] = fileparts(imds.Files{i});
    new_class_matrix(i,1) = strsplit(img_name,fsp); 
end

%Get image classifications
new_class_matrix(:,2) = cellstr(imds.Labels(:)); 

%Sort imds by image name
sort_help(:,1) = strtok(new_class_matrix(:,1),'_'); sort_help = str2double(sort_help);
new_class_matrix(:,3)=num2cell(sort_help);
new_class_matrix = sortrows(new_class_matrix,3);

%% Adjust & Write new_XLSX
new_matrix_helper = readtable(new_XLSX,'Sheet','Sheet1');
new_matrix_helper = table2cell(new_matrix_helper(1:end-4,2:end));

%Adjust matrix with new classifications
new_class_matrix = reshape(new_class_matrix(:,2), size(new_matrix_helper,2), size(new_matrix_helper,1))';

%Write changes to new file
xlswrite(new_XLSX,new_class_matrix,'Sheet1','B2');

%Delete old file
delete(old_XLSX);

clearvars -except fsp XLSXdir;
end
catch err
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%% ERROR while creating new.xlsx %%%%')
        disp(err.message)
        disp('Remove Adjusted_ from .xlsx filename')
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        return
end

t1 = toc;   
disp(['%% batch_ReClassifier: Finished everything successfully in ' num2str(t1) ' seconds %%']) 