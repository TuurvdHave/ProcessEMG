function imgCreator(parLoop,header,body,body2,EMG_amps,dirImages)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    imgCreator for the EMG Classifier                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <Riad.Akhundov@uon.edu.au>
%_____________________________________________________________________________
%Function for creating EMG images based on number of c3d files in acquisition.

%% Check which loop to run
fsp = filesep;
if parLoop == 1
    %Parallel loop for more than two .c3d files
    parfor i=1:length(header)
        if strcmp(header(3,i),'No EMG detected') == 1
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['%%%% ERROR in ' header{2,i} ' %%%%'])
            disp('%%%%      No EMG detected in this c3d %%%%')
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        else
            f = figure('visible', 'off');
            plot(body(:,i),'color',[0, 0.4470, 0.7410]); %Bandpassed signal
            hold on;
            plot(body2(:,i),'m','LineWidth',4); %Envelope

            hline = line([0,sum(~isnan(body(:,i)))],[0,0]); %Reference line at 0
            hline.Color = 'm';
            hline.LineWidth = 2;

            axis tight;
            ylim([-EMG_amps(1,i),EMG_amps(1,i)]); %Limit amplitude here
            set(gca,'xtick',[],'ytick',[]);
            set(gca, 'Box', 'off');
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.27 2.27]);
            set(gca,'LooseInset',get(gca,'TightInset'));
            set(gca,'YColor','none','XColor','none')

            %Creates images all over in the dir even if they already exist 
            print (strcat(dirImages, fsp, header(4,i), '.jpg'), '-djpeg', '-r100'); 
            close(f); 
        end
    end

else
    
    %Normal loop for less than three .c3d files
    for i=1:length(header)
        if strcmp(header(3,i),'No EMG detected') == 1
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            disp(['%%%% ERROR in ' header{2,i} ' %%%%'])
            disp('%%%%      No EMG detected in this c3d %%%%')
            disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        else
            f = figure('visible', 'off');
            plot(body(:,i),'color',[0, 0.4470, 0.7410]); %Bandpassed signal
            hold on;
            plot(body2(:,i),'m','LineWidth',4); %Envelope

            hline = line([0,sum(~isnan(body(:,i)))],[0,0]); %Reference line at 0
            hline.Color = 'm';
            hline.LineWidth = 2;

            axis tight;
            ylim([-EMG_amps(1,i),EMG_amps(1,i)]); %Limit amplitude here
            set(gca,'xtick',[],'ytick',[]);
            set(gca, 'Box', 'off');
            set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.27 2.27]);
            set(gca,'LooseInset',get(gca,'TightInset'));
            set(gca,'YColor','none','XColor','none')

            %Creates images all over in the dir even if they already exist 
            print (strcat(dirImages, fsp, header(4,i), '.jpg'), '-djpeg', '-r100'); 
            close(f); 
        end
    end
end
end
