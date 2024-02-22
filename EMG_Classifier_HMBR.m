%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              EMG Classifier                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <Riad.Akhundov@uon.edu.au>

% Adapted by Miel Willems and van der Have Tuur to fit the Human movement
% research group Leuven

%Main script for classifying EMG quality: Read a participant-folder containing
%c3d files, classify the EMG quality, save results in xlsx files and sort 
%created EMG images according to their class.

%%%Disclaimer:
% By downloading, accessing or using this software, you signify your assent to this disclaimer. 
% The contents of this app, including without limitation, all data, information, text, links and other
% materials are meant to be used for informational purposes only. We do not take responsibility for 
% decisions taken by the user based solely on the information provided in this software.
% Nothing contained in this software is intended to create a physician-patient relationship, to replace
% the servicesof a licensed, trained physician or health professional or to be a substitute for medical 
% advice of a licensed physician or trained health professional. You should consult a licensed physician 
% in all matters relating to your health. 
% You hereby agree that you shall not make any health or medical related decision based in whole or in 
% part on anything contained in this software.

%%%Requirements: 
%1) MATLAB 2018b or newer (newest MATLAB version is strongly recommended)
%2) Neural Network Toolbox (renamed to Deep Learning Toolbox since 2018a)
%3) Parallel Computing Toolbox (also called Distributed Computing Toolbox)

%Version: v0.99a
clc; close all; clearvars -except EMG_Classifier_AlexNet_v4d;

%Select folder with .c3d Files
fsp = filesep;
folder = uigetdir(path,'Select Participant Folder Containing .c3d Files');
fileList = dir([folder, fsp,'**',fsp,'*.c3d']);
[~,participant] = fileparts(folder);

%Check if selected folder contains .c3d Files
if isempty(fileList) == 1
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%% ERROR: No .c3d files in selected folder %%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    return
end

%Asking input of the user 

Firstanswer = inputdlg({'Low pass filter (Hz)','Band pass filter (Hz)','Manual crop?','White EMG used?','Blue EMG used?','Outside?','Filter out ECG?'},'Analyses',[1 35],{'6','[20 400]','Yes/No','Not used/Used','Not used/Used','Yes/No','Yes/No'});

if strcmpi(Firstanswer{4,1},'used') 
submuscles = inputdlg({'Muscle channel White 1?','Muscle channel White 2?','Muscle channel White 3?',...
    'Muscle channel White 4?','Muscle channel White 5?','Muscle channel White 6?','Muscle channel White 7?','Muscle channel White 8?','Muscle channel White 9?',...
    'Muscle channel White 10?','Muscle channel White 11?','Muscle channel White 12?','Muscle channel White 13?','Muscle channel White 14?','Muscle channel White 15?',...
    'Muscle channel White 16?'},'Analyses',[1 35],{'Not used','Not used','Not used','Not used',...
    'Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used'});
else 
    submuscles = {'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used'};
end 
if strcmpi(Firstanswer{5,1},'used') 
    submuscles(17:32) = inputdlg({'Muscle channel Blue 1?','Muscle channel Blue 2?','Muscle channel Blue 3?',...
    'Muscle channel Blue 4?','Muscle channel Blue 5?','Muscle channel Blue 6?','Muscle channel Blue 7?','Muscle channel Blue 8?','Muscle channel Blue 9?',...
    'Muscle channel Blue 10?','Muscle channel Blue 11?','Muscle channel Blue 12?','Muscle channel Blue 13?','Muscle channel Blue 14?','Muscle channel Blue 15?',...
    'Muscle channel Blue 16?'},'Analyses',[1 35],{'Not used','Not used','Not used','Not used',...
    'Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used','Not used'});
else 
    submuscles(17:32) = {'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used';'Not used'};
end 
a = 1; 
for i = 1:32
    if strcmpi(submuscles(i),'not used')
    else 
    channels_used(a) = i; 
    muscles(a) = submuscles(i);
    a = a +1;
    end 
end 

if strcmpi(Firstanswer{7,1},'yes')
   ecg_answer = inputdlg(muscles,'ECG-filtering per muscle, answer with yes or no',[1 35]);
end 

%% Getting Data & Requirements Checks
%Change the current folder to the folder of this m-file (if user does "Add to Path" instead of "Change Folder")
if(~isdeployed) 
  cd(fileparts(which(mfilename)));
end

%Check for requiured MATLAB version & Toolboxes
toolbox_helper(1:5,1) = {'Correct MATLAB Version';'License - Neural Network';'License - Parallel Computing';...
    'Installed - Neural Network';'Installed - Parallel Computing'};
toolbox_helper(1:5,2) = num2cell([~verLessThan('matlab','9.5'); license('test','Neural_Network_Toolbox');...
    license('test','Distrib_Computing_Toolbox');~isempty(ver('nnet'));~isempty(ver('distcomp'))]);

if any(cell2mat(toolbox_helper(:,2)) == 0)
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('%%%% ERROR: MATLAB version 2018b or newer required %%%%')
    disp('%%%%        Neural Network Toolbox required        %%%%')
    disp('%%%%        Parallel Computing Toolbox required    %%%%')
    disp('%%%%                                               %%%%')
    disp('%%%%       - Check toolbox_helper in Workspace -   %%%%')
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    return
end

%Load trained AlexNet
if ~exist('EMG_Classifier_AlexNet_v4d','var') == 1
    str = mfilename('fullpath');
    cd(str(1:end-19))
    addpath([pwd, fsp,'bin'])
    addpath([pwd, fsp,'ecg-removal'])
    load ([pwd, fsp,'bin',fsp,'EMG_Classifier_AlexNet_v4d']) 
end
tic;
    addpath([pwd, fsp,'bin'])
    addpath([pwd, fsp,'bin/ecg-removal'])
    folder2 = fileparts(which([pwd, fsp,'bin/ecg-removal/ecg_removal.m']));
    addpath(genpath(folder2));
%Turn unnecessary warnings off 
warning('off', 'btk:ReadAcquisition')
warning('off', 'MATLAB:MKDIR:DirectoryExists')
warning('off', 'MATLAB:xlswrite:AddSheet')
warning('off', 'MATLAB:rmpath:DirNotFound')
warning('off', 'MATLAB:DELETE:FileNotFound')
warning('off', 'MATLAB:hg:AutoSoftwareOpenGL')

%% Directories Creation
%Create directories in previously chosen folder
dirClassifier = [folder, fsp, 'Classifier_', participant];
mkdir(dirClassifier);

dirImages = [dirClassifier, fsp, 'Images'];
if exist(dirImages,'dir')
    A = dir([dirImages, fsp, '1_Good']);
        for k = 1:length(A)
            delete([dirImages, fsp, '1_Good',  '\' A(k).name])
        end
    clearvars A 
    A = dir([dirImages, fsp, '2_Usable']);
        for k = 1:length(A)
            delete([dirImages, fsp, '2_Usable',  '\' A(k).name])
        end
    clearvars A 
    A = dir([dirImages, fsp, '3_Noise']);
        for k = 1:length(A)
            delete([dirImages, fsp, '3_Noise',  '\' A(k).name])
        end
    clearvars A 
    A = dir([dirImages, fsp, '4_NoSignal']);
        for k = 1:length(A)
            delete([dirImages, fsp, '4_NoSignal',  '\' A(k).name])
        end
    clearvars A 
else
mkdir(dirImages);
mkdir([dirImages, fsp, '1_Good']);
mkdir([dirImages, fsp, '2_Usable']);
mkdir([dirImages, fsp, '3_Noise']);
mkdir([dirImages, fsp, '4_NoSignal']);
end 

dirXLSX = [dirClassifier, fsp, 'XLSX'];
mkdir(dirXLSX);
addpath([pwd, fsp,'bin']); %Add path to required functions
addpath([pwd, fsp,'bin', fsp,'btk']); %Add path to BTK

%XLSX directories
old_XLSX = [pwd, fsp, '_Custom Scripts', fsp, 'Classifications_Template.xlsx'];
new_XLSX = [dirXLSX, fsp,'Classifications_', participant, '.xlsx'];
delete(new_XLSX); %Delete old xlsx version to overwrite it
copyfile (old_XLSX,new_XLSX);

%% Parallel Computing Check and Diary
%Get rid of static trials
for s=length(fileList):-1:1  
    if contains(upper(fileList(s).name),'STATIC') %Static trials should be named "static..."
        fileList(s) = [];
    elseif contains(upper(fileList(s).name),'CAL') %Static trials should be named "cal..."
        fileList(s) = [];
    end
end

numTrials = length(fileList); 
numEMG(1:numTrials,1) = NaN; %Preallocation

%Parallel computing check
if numTrials>2
    %Start parallel computing pool
    poolobj = gcp('nocreate'); 
    if isempty(poolobj)
        poolobj = parpool;
    end
    parLoop = 1; %Parallel check variable
else
    parLoop = 0;
end

%Create diary 
diaryName = ['log_', participant];
diary off; fclose( 'all' ); delete(diaryName);  %Against movefile errors when aborting and restarting the script
diary(diaryName); 
diary on
disp('%% Started reading all c3d %%');disp('%')

%% Take & Filter Data from .c3d Files
for s=1:numTrials
    try
    %Get Analog Data into mat
    h=btkReadAcquisition(fullfile(fileList(s).folder, fileList(s).name));
    [matFiles(s).analogs, matFiles(s).analogsInfo] = btkGetAnalogs(h);
    btkDeleteAcquisition(h);
    c3dName(s,1) = cellstr(fileList(s).name(1:end-4)); %Grab c3d names for later

    %Turn structs into cells
    matFiles(s).analogsInfo.label2 = struct2cell(matFiles(s).analogsInfo.label);
    matFiles(s).analogsInfo.units2 = struct2cell(matFiles(s).analogsInfo.units);
    matFiles(s).EMGs = struct2cell(matFiles(s).analogs);

    %Find where EMGs are
    label2 = ~(contains(upper(matFiles(s).analogsInfo.label2),'FORCE')+...  %Assuming EMG is not called "Channel"
        contains(upper(matFiles(s).analogsInfo.label2),'MOMENT')+...        %Assuming EMG is not called "Timing_gate"
        contains(upper(matFiles(s).analogsInfo.label2),'RAW')+...           %Assuming EMG is not called "TG" for Timing_gate
        contains(upper(matFiles(s).analogsInfo.label2),'ACCEL'));           %Assuming EMG is not called "Sync"
                                                                            %<-- Add further exclusions here if required
    
    units2 = contains(upper(matFiles(s).analogsInfo.units2),'V');           %EMG units are always in V 
    locEMG = find(units2==1 & label2==1);                                   %Find EMG locations
    locEMG = locEMG';
    if strcmpi(Firstanswer{6,1},'No')
    locEMG(:,33:34) = [];
    end 
    %select only channels that were used in acquisition
    for channel = 1:size(channels_used,2)
         locEMG_upd(:,channel) = locEMG(:,channels_used(channel));   
    end
    
    %Combine important info in struct
    matFiles(s).EMGs = matFiles(s).EMGs(locEMG_upd)';
    matFiles(s).RawEMGs = horzcat(matFiles(s).EMGs{:});

    % Cut out the gaitcycle manually
    if strcmpi((Firstanswer{3,1}),'yes')
    
    timelist = (1:length( matFiles(s).RawEMGs))';
    timelist = timelist./1000;
    Secondanswer = inputdlg({'start crop time in sec','end crop time in sec'},char(c3dName{s,1}),[1 35],{'0','4.5'});
    start_ind = find(timelist(:,1)==str2num(Secondanswer{1,1}));
    stop_ind = find(timelist(:,1)==str2num(Secondanswer{2,1}));

    matFiles(s).RawEMGs = matFiles(s).RawEMGs(start_ind:stop_ind,:);
    else %do nothing
    end 
    
    matFiles(s).RawEMGs = matFiles(s).RawEMGs./10^-3; %scaling factor 
    matFiles(s).EMGlabels = matFiles(s).analogsInfo.label2(locEMG_upd)';
    matFiles(s).EMGlabels = muscles;
    matFiles(s).EMGs = vertcat(matFiles(s).EMGlabels,matFiles(s).EMGs);
    matFiles(s).EMGs{1,end+1} = [fileList(s).folder, fsp, fileList(s).name];
    matFiles(s).EMGs{2,end} = matFiles(s).analogsInfo.frequency;

    %Filtering
    Rate = matFiles(s).EMGs{2,end};
    [b,a] = butter(4,str2num(Firstanswer{2,1})./(Rate/2),'bandpass');
    [b2,a2] = butter(4,str2num(Firstanswer{1,1})./(Rate/2),'low');

    for col=1:size(matFiles(s).RawEMGs,2)
        matFiles(s).RawEMGs(:,col) = matFiles(s).RawEMGs(:,col) - mean(matFiles(s).RawEMGs(:,col)); %Raw offset-corrected data
        if strcmpi(Firstanswer{7,1},'yes') && strcmpi(ecg_answer{col,1},'yes')
            matFiles(s).RawEMGs(:,col) = ecg_removal(matFiles(s).RawEMGs(:,col)'); % ECG-removal
        end 
        matFiles(s).FiltEMGs(:,col) = abs(filtfilt(b,a,matFiles(s).RawEMGs(:,col))); %Bandpass filtered data
        matFiles(s).Envelopes(:,col) = filtfilt(b2,a2,(matFiles(s).FiltEMGs(:,col))); %Lowpass filtered data
    end
    
    numEMG(s,1) = length(matFiles(s).EMGlabels); %Get number of EMGs ix`n this c3d
    catch err
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(['%%%% ERROR in ' fileList(s).name ' %%%%'])
        disp(err.message)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
    EMG.(fileList(s).name(1:end-4)) = matFiles(s);
end
[numEMGmax, locEMGmax] = max(numEMG); %Find maximum number of EMG for later sorting

userview = memory; mem1 = userview.MemUsedMATLAB/1024^3; %Memory usage
t1 = toc;
disp(['%% Finished reading all c3d in ' num2str(t1) 'sec, using ' num2str(mem1) 'GB of RAM' ' %%']);disp('%')

%% Data Matrix Creation for XLSX
disp('%% Started data manipulation %%');disp('%')

%Preallocations
headerC = cell(4,numTrials);
bodyC = cell(1,numTrials);
body2C = cell(1,numTrials);
bodyC_pl = cell(1,numTrials);
body2C_pl = cell(1,numTrials);
EMG_amps = zeros(numTrials,numEMGmax);

%Check for EMG consistency between trials
if numel(unique(numEMG))==1
    numEMGcheck=1;
else
    numEMGcheck=2;
end

%Create matrix 
for s=1:numTrials
    headerC{1,s} = strings(1,1);
    headerC{2,s} = strings(1,1);
    headerC{3,s} = strings(1,1);
    headerC{4,s} = strings(1,1);
        
    if isempty(matFiles(s).EMGlabels)  
            headerC{1,s}(1,1) = [fileList(s).folder, fsp];  %If no EMG detected
            headerC{2,s}(1,1) = [fileList(s).name]; 
            headerC{3,s}(1,1) = 'No EMG detected'; 
    elseif numEMGcheck ==1                                  %If some EMG detected and all trials have same number of EMGs
        for i=1:length(matFiles(s).EMGlabels)                
            headerC{1,s}(1,i) = [fileList(s).folder, fsp];
            headerC{2,s}(1,i) = [fileList(s).name];
            headerC{3,s}(1,i) = matFiles(s).EMGlabels{1,i};
        end                     
            bodyC{1,s} = matFiles(s).FiltEMGs;
            body2C{1,s} = matFiles(s).Envelopes;
    else                                                    %If trials have different number of EMGs
            bodyC_pl{1,s} = matFiles(s).FiltEMGs;           %Placeholder for sorting out missing EMGs
            body2C_pl{1,s} = matFiles(s).Envelopes;
            
            bodyC{1,s} = zeros([size(matFiles(s).FiltEMGs,1), numEMGmax],'double');
            body2C{1,s} = bodyC{1,s};
            
        for i=1:numEMGmax
            headerC{1,s}(1,i) = [fileList(s).folder, fsp];
            headerC{2,s}(1,i) = [fileList(s).name];
            headerC{3,s}(1,i) = matFiles(locEMGmax).EMGlabels{1,i}; %Labels without missing EMGs   

        end
            
        for i=1:length(matFiles(s).EMGlabels)                
            headerC{4,s}(1,i) = matFiles(s).EMGlabels{1,i}; %Labels with missing EMGs 
            bodyC{1,s}(:,find(ismember(headerC{3,s},headerC{4,s}{1,i}))) = bodyC_pl{1,s}(:,i); %Important line: Sorting EMGs
            body2C{1,s}(:,find(ismember(headerC{3,s},headerC{4,s}{1,i}))) = body2C_pl{1,s}(:,i); %Ignore warning
        end     
    end
end

%Adjust header
header(1,:) = horzcat(headerC{1,1:end});
header(2,:) = horzcat(headerC{2,1:end});
header(3,:) = horzcat(headerC{3,1:end});
header(4,:) = strings(1,1);

for i=1:length(header)
    if strcmp(header(3,i),'No EMG detected') == 1
        header(4,i) = 'No images created';
    else
        imgName = strcat(num2str(i), '_', header(3,i), '_', header(2,i));
        header(4,i) = imgName{1,1}(1:end-4);
    end
end

header_cell = cellstr(header);

%Adjust body
LengthCells=cellfun('size',bodyC,1);
maxLengthCell=max(cellfun('size',bodyC,1)); %Finding the longest vector in the cell array

for i=1:length(bodyC)
    
    for j=cellfun('size',bodyC(i),1)+1:maxLengthCell
        if isempty(matFiles(i).EMGlabels)
            bodyC{i}(1:maxLengthCell,1)=NaN;
            body2C{i}(1:maxLengthCell,1)=NaN;
        else
            bodyC{i}(j,1:numEMGmax)=NaN;   %NaN-pad the elements in each cell array with a length shorter than the maxlength
            body2C{i}(j,1:numEMGmax)=NaN;          
        end
    end
   
end
body=cell2mat(bodyC); %Combining FiltEMGs into double matrix
body2=cell2mat(body2C); %Combining Envelopes into double matrix

userview = memory; mem2 = userview.MemUsedMATLAB/1024^3;
t2 = toc;
disp(['%% Finished data manipulation in ' num2str(t2) 'sec, using ' num2str(mem2) 'GB of RAM' ' %%']);disp('%')

% Make Images
disp('%% Started creating images %%');disp('%')
%Calculate Amplitudes of FiltEMGs for Image-creation
for i=1:length(matFiles)
    EMG_amps(i,:) = max(rmoutliers(cell2mat(bodyC(1,i)),'movmean',(length(cell2mat(bodyC(1,i)))/5)));
end
EMG_amps(EMG_amps == 0) = NaN; EMG_amps = mean(EMG_amps,1,'omitnan'); EMG_amps = repmat(EMG_amps,1,numTrials);   

%Create Images
imgCreator(parLoop,header,body,body2,EMG_amps,dirImages)

userview = memory; mem3 = userview.MemUsedMATLAB/1024^3;
t3 = toc;
disp(['%% Finished creating images in ' num2str(t3) 'sec, using ' num2str(mem3) 'GB of RAM' ' %%']);disp('%')

%% Classification
disp('%% Started classifying %%');disp('%')

%Adjust header#2
try
header_cell2 = header_cell'; %Turn to cell to combine later
indices = find(strcmp(header(3,:),'No EMG detected') == 1); %Find missing EMG
header_cell2(indices,:) = []; %Delete missing EMG

%Setup the imageDatastore
imds = imageDatastore(dirImages, 'IncludeSubfolders', false);
new_class_matrix = cell(length(imds.Files),7);

%Classify images
[labels, percentages] = classify(EMG_Classifier_AlexNet_v4d, imds); %The classification line
percentages_round = floor(percentages*10000)/100;

%Get image names
for i=1:length(imds.Files)
    [~,img_name] = fileparts(imds.Files{i});
    new_class_matrix(i,1) = strsplit(img_name,fsp); 
end

%Get image classifications
new_class_matrix(:,2) = cellstr(labels); 
new_class_matrix(:,3:end-1) = num2cell(percentages_round);

%Sort imds by image name
sort_help(:,1) = strtok(new_class_matrix(:,1),'_'); sort_help = str2double(sort_help);
new_class_matrix(:,end) = num2cell(sort_help);
new_class_matrix = sortrows(new_class_matrix,7);

%Put results in one big cell
labels_cell = new_class_matrix(:,2);
percentages_cell = num2cell(percentages_round); 
class_matrix = horzcat(header_cell2,new_class_matrix(:,2:6)); %Cell with header, class & rounded percentages
catch err
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%% ERROR while classifying %%%%')
        disp(err.message)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

for i = 1 : size(fileList,1)
    ch = size(channels_used,2);
    da = class_matrix((i*ch-ch+1):(ch*i),5);
    for ii = 1: size(da,1)
    EMG.(fileList(i).name(1:end-4)).Class(ii) = str2double(da{ii,1}(1));
    end 
end 

%%  save .mat and xls 
% mat 
save(fullfile(fileList(s).folder,'EMG.mat'),'EMG')
% raw EMG in xls 
tr = fields(EMG);
for i = 1:size(fields(EMG),1)
    table(EMG.(tr{i}).RawEMGs);
    writetable(ans,fullfile(fileList(s).folder,'Raw.xlsx'),'Sheet',tr{i});
end 
table(EMG.(tr{i}).EMGlabels);
writetable(ans,fullfile(fileList(s).folder,'Raw.xlsx'),'Sheet','EMGlables');
% envelopes EMG in xls 
for i = 1:size(fields(EMG),1)
    table(EMG.(tr{i}).Envelopes);
    writetable(ans,fullfile(fileList(s).folder,'Envelopes.xlsx'),'Sheet',tr{i});
end 
table(EMG.(tr{i}).EMGlabels);
writetable(ans,fullfile(fileList(s).folder,'Envelopes.xlsx'),'Sheet','EMGlables');


userview = memory; mem4 = userview.MemUsedMATLAB/1024^3;
t4 = toc;
disp(['%% Finished classifying in ' num2str(t4) 'sec, using ' num2str(mem4) 'GB of RAM' ' %%']);disp('%')

%% Saving Results in XLSX:
disp('%% Saving results %%');disp('%')

try
%Making the XLSX varnames    
lenEMG = length(numEMG);
lenEMG_nonzero = sum(numEMG>0);

class_matrix_xlsx = cell(lenEMG+5,numEMGmax+1);

separator = cellstr('--------------');
noEMG = cellstr('-No EMG-');

class_matrix_xlsx(:) = noEMG;
class_matrix_xlsx{1,1} = 'c3d/Muscle';
class_matrix_xlsx(2:lenEMG+1,1) = c3dName;

class_matrix_xlsx(lenEMG+2,1) = cellstr('--Mode--');
class_matrix_xlsx(lenEMG+2,2:end) = separator;
class_matrix_xlsx{lenEMG+3,1} = 'EMG-driven';
class_matrix_xlsx{lenEMG+4,1} = 'EMG-hybrid';
class_matrix_xlsx{lenEMG+5,1} = 'EMG-assisted';

class_matrix_xlsx(1,2:numEMGmax+1) = matFiles(locEMGmax).EMGlabels; %Muscle names  

%Filling in the classifications
k=1;
for i = 1:lenEMG 
    if numEMG(i) == 0
        class_matrix_xlsx(i+1,2:numEMGmax+1) = noEMG;
    else
        class_matrix_xlsx(i+1,2:numEMGmax+1) = class_matrix(k:numEMGmax+k-1,5)';
        k = k+numEMGmax;
    end
        
end

%Filling in the correct mode
for j=2:numEMGmax+1
    numNoise = sum(ismember(class_matrix_xlsx(2:lenEMG+1,j), {'3_Noise', '4_NoSignal'}));
    numUsable = sum(ismember(class_matrix_xlsx(2:lenEMG+1,j), {'2_Usable'}));
    if numNoise == 0 && numUsable == 0
        class_matrix_xlsx{lenEMG+3,j} = 'all good';
        class_matrix_xlsx{lenEMG+4,j} = 'all good';
        class_matrix_xlsx{lenEMG+5,j} = 'all good';
       
    elseif numNoise > floor(lenEMG_nonzero*0.5)
        class_matrix_xlsx{lenEMG+3,j} = 'remap';
        class_matrix_xlsx{lenEMG+4,j} = 'synthesize';
        class_matrix_xlsx{lenEMG+5,j} = 'synthesize';
        
    else
        class_matrix_xlsx{lenEMG+3,j} = 'remap';
        class_matrix_xlsx{lenEMG+4,j} = 'remap';
        class_matrix_xlsx{lenEMG+5,j} = 'adjust';
    end
end

%Write XLSX
xlswrite(new_XLSX,class_matrix_xlsx,'Sheet1');
xlswrite(new_XLSX,class_matrix,'Sheet2');
catch err
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%% ERROR while creating .xlsx %%%%')
        disp(err.message)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end


%% Saving Results by Sorting Images
try
for i=1:length(labels_cell)
    if strcmp(labels_cell{i,1},'1_Good') == 1
        movefile([dirImages, fsp, header_cell2{i,4}, '.jpg'], [dirImages, fsp, '1_Good']);  
        
    elseif strcmp(labels_cell{i,1},'2_Usable') == 1
        movefile([dirImages, fsp, header_cell2{i,4}, '.jpg'], [dirImages, fsp, '2_Usable']);
        
    elseif strcmp(labels_cell{i,1},'3_Noise') == 1
        movefile([dirImages, fsp, header_cell2{i,4}, '.jpg'], [dirImages, fsp, '3_Noise']);
        
    else
        movefile([dirImages, fsp, header_cell2{i,4}, '.jpg'], [dirImages, fsp, '4_NoSignal']);
    end
end
catch err
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%% ERROR while sorting images %%%%')
        disp(err.message)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end

%% Close and Move Diary
userview = memory; mem5 = userview.MemUsedMATLAB/1024^3;
t5 = toc;
disp(['%% EMG_Classifier: Finished everything successfully in ' num2str(floor(t5/60)) ' minutes and ' num2str(rem(t5,60))...
    ' seconds, using ' num2str(mem5) 'GB of RAM' ' %%'])

diary off
try
    movefile(diaryName, [dirClassifier, fsp, diaryName, '.txt']); %Creates diary in folder selected at the beginning
catch err
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp('%%%% ERROR while moving diary %%%%')
        disp(err.message)
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
end
