%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              EMG Normalizer                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code will normalize your EMG to the maximum values found in the
% selected c3d-files. This could be meeasured MVC's but also other motions.
% It is recommended to only choose files that are classfied usable or good
% when normalizing as peaks in these data could mess up the normalization. 

clear all; close all 

%Select folder with .c3d Files
[file, path]= uigetfile('*.mat',  'Select the mat structure containing the EMG signals . . .',  'MultiSelect', 'off');
load(fullfile(path,file));
[c3dFilenames,path] = uigetfile('*.c3d',  'Select c3d files you want to normalize to . . .',  'MultiSelect', 'on');

% creating matrix containing maximum values found for the C3D files
% selected
lastMVC = [];
for i = 1 : size(c3dFilenames,2)
    MVC = max(EMG.(c3dFilenames{i}(1:end-4)).Envelopes,[],1);
    if sum(size(lastMVC)) == 0
        lastMVC = MVC;
    else
        for tt = 1 : length(lastMVC)
            if lastMVC(tt) < (MVC(tt))
                lastMVC(tt) = (MVC(tt));
            end
        end
    end
end

% normalizing to the maximum value defined in lastMVC
trials = fields(EMG);
for i = 1:size(fields(EMG),1)
    for c = 1:size(EMG.(trials{i}).Envelopes,2)
        EMG.(trials{i}).Norm(:,c) = EMG.(trials{i}).Envelopes(:,c)./lastMVC(c);
    end 
end 

