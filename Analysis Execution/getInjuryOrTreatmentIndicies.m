function [InjuryIndicies,AnalysisRegions, plotProps] = getInjuryOrTreatmentIndicies(RecordingMetrics,plxFileName,sw,plotProps,ProcessedData)
%This function allows the user to specify where treatments were
%administered during the recording and plots recording outputs independent
%of correlograms with the labeled treatments and analysis regions

%% First note where in the recording an injury or treatment occurred (or whether there was no injury or treatment)
if ~exist([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'file') || ~exist([pwd '\Output Files\' plxFileName '\InjuryIndicies.mat'],'file') || sw.setAnalysisRegions == 1 || sw.setInjuryRegions == 1 %if you have not yet processed this data or would like to re-process the data...

    if ~exist([pwd '\Output Files\' plxFileName],'dir'); mkdir([pwd '\Output Files\' plxFileName]); end %make the folder to store matlab outputs for the given recording if it does not already exist

    if sw.setInjuryRegions == 1 || ~exist([pwd '\Output Files\' plxFileName '\InjuryIndicies.mat'],'file') %If want or need to label treatments/injuries
        figure() %make a plot so that you can identify where the injury or treatment occurred
        h1=plot(RecordingMetrics.minBinsSPM,RecordingMetrics.meanSPM);
        xlabel('Time (min)')
        ylabel('Events per Minute')
        set(gca,'FontSize',plotProps.FigureFontSize)
        xlim([RecordingMetrics.minBinsSPM(1) RecordingMetrics.minBinsSPM(end)])

        disp(' ')
        %Based on the plot, you can usually tell when an injury or treatment was delivered.  Get the start and end time in minutes
        multipleInjIndicator = input('Are there multiple injuries and/or treatments in this recording?  Type 1 if yes, 0 otherwise (including no treatment/injury). ');

        if multipleInjIndicator ~= 1 %if only one or fewer injuries or treatments
            InjuryStartMin = input('     Please type the time (x value of the plot) an injury/treatment started being administered.  If no injury/treatment, type 0.  If injury/treatment was administered before the start of the recording, type -1.  ');
            if InjuryStartMin>1 %if the injury was administered during the recording
                InjuryEndMin = input('     Please type the time (x value of the plot) the injury/treatment finished being administered. ');
            else %otherwise, the injury/treatment was given before the recording started, so injury/treatment end is also before the recording started
                InjuryEndMin = -1;
            end
            InjuryLabel = input('     Please type a label of what the injury or treatment was (i.e. pH shock, bicuculline, etc...).  Enter 0 if none. ','s');
            InjuryLabels = {}; InjuryLabels = [InjuryLabels, InjuryLabel, 'q']; %add an extra blank entry ('q') in the label list so matlab can index labels correctly

        else %otherwise, multiple injuries or treatments were administered

            InjuryStartMin = []; InjuryEndMin = []; %initialize injury/treatment time arrays
            InjuryLabels = { }; %initialize injury/treatment label array
            injuriesDone = 0; %initialize an indicator that determines whether or not to keep selecting injuries

            while injuriesDone == 0  %while still selecting injuries/treatments...
                InjuryStartMina = input('     Please type the time (x value of the plot) an injury/treatment started being administered.  If injury/treatment was administered before the start of the recording, type -1.  ');
                if InjuryStartMina>1 %if the injury was administered during the recording
                    InjuryEndMina = input('     Please type the time (x value of the plot) the injury/treatment finished being administered. ');
                else %otherwise, the injury/treatment was given before the recording started, so injury end was also before the recording started
                    InjuryEndMina = -1;
                end
                InjuryLabel = input('     Please type a label of what the injury or treatment was (i.e. pH shock, bicuculline, etc...). ','s');

                hold on
                injStartIdx = max(1,find(RecordingMetrics.minBinsSPM<max(2,InjuryStartMina),1,'last'));
                injEndIdx =  max(1,find(RecordingMetrics.minBinsSPM<max(2,InjuryEndMina),1,'last'));
                h2=plot(RecordingMetrics.minBinsSPM(injStartIdx:injEndIdx),RecordingMetrics.meanSPM(injStartIdx:injEndIdx),'m','LineWidth',2);
                hold off
                legend([h1 h2],'Entire Recording','Selected Injuries/Treatments')
                xlim([RecordingMetrics.minBinsSPM(1) RecordingMetrics.minBinsSPM(end)])

                injuryMistake = input('     Was there a mistake in any of the entries for this treatment/injury?  Type 1 if yes, 0 if no. If yes, do not fret.  You can re-enter the correct information for this injury/treatment. ');
                if injuryMistake == 0 %if there are no mistakes, save the injury/treatment
                    InjuryStartMin = [InjuryStartMin; InjuryStartMina];
                    InjuryEndMin = [InjuryEndMin; InjuryEndMina];
                    InjuryLabels = [InjuryLabels, InjuryLabel];

                    injuriesDone = input('Are you finished selecting injuries/treatments?  Type 1 if yes, 0 if no. ');

                    %If there was a mistake, the code will just repeat the loop again without saving the incorrect information
                end

            end %end pick injuries/treatments

        end %end check for multiple injuries
        close(gcf)

        %Sort injuries/treatments chronologically
        [InjuryStartMin, sortInds] = sortrows(InjuryStartMin,1); %sort the region entries so first entry = earliest in the recording
        InjuryLabels = InjuryLabels(sortInds); %sort the region labels too
        InjuryEndMin = InjuryEndMin(sortInds); %sort the injury ends too

        %The units of the raster are in seconds, convert injury minutes to seconds
        InjuryStartSec = InjuryStartMin.*60;
        InjuryEndSec = InjuryEndMin.*60;

        %Save the injury start and end in the InjuryIndicies structure
        InjuryIndicies.InjuryStartMin = InjuryStartMin;
        InjuryIndicies.InjuryEndMin = InjuryEndMin ;
        InjuryIndicies.InjuryStartSec = InjuryStartSec;
        InjuryIndicies.InjuryEndSec =  InjuryEndSec;
        InjuryIndicies.InjuryLabels = InjuryLabels;
        if InjuryStartMin~=0 %if there was an injury/treatment
            InjuryIndicies.NumberOfInjuriesOrTreatments = size(InjuryEndMin,1);
        else %otherwise no injury/treatment
            InjuryIndicies.NumberOfInjuriesOrTreatments = 0;
        end

        save([pwd '\Output Files\' plxFileName '\InjuryIndicies.mat'],'InjuryIndicies')

    else %otherwise load pre-selected injury labels and times

        load([pwd '\Output Files\' plxFileName '\InjuryIndicies.mat'],'InjuryIndicies');

    end %end check whether to set or load injury/treatment labels


    %% Now select regions of the data to analyze
    if ~exist([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'file') || sw.setAnalysisRegions == 1
        
        if exist([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'file') && sw.LabelJunkRegions == 0 && sw.LabelBadSignals == 0 %create the folder to save code outputs if this has not already been done
            disp(' ')
            AddToRegionsIndicator = input('Do you wish to manually add to previously-selected data regions?  If yes, type 1.  If no, type 0 and all regions will be selected from scratch.  ');
        else
            AddToRegionsIndicator = 0;
        end

        if ~exist([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'file') || AddToRegionsIndicator == 0

            %Once the injury/treatment is identified, select regions of the data to analyze
            disp(' ')
            DivideIndicator = input(['Do you wish to select data regions to analyze by hand?  If yes, type 1.  If no, type 0 and the algorithm will make ' num2str(((plotProps.automatedRegionBinLength-rem( max(RecordingMetrics.minBinsSPM),plotProps.automatedRegionBinLength))+ max(RecordingMetrics.minBinsSPM))/plotProps.automatedRegionBinLength) ', ' num2str(plotProps.automatedRegionBinLength) ' min regions from recording start to finish.  ']);
            if DivideIndicator == 1 %if you want to analyze several sections of the recording seperately...

                doneDividing = 0; %make a switch that will set to 1 once you have finished selecting regions of the data to analyze

                figure()
                hold on
                h1=plot(RecordingMetrics.minBinsSPM,RecordingMetrics.meanSPM,'k');
                hold off
                xlabel('Time (min)')
                ylabel('Events per Minute After Injury or Treatment')
                xlabel('Time (min)')
                ylabel('Events per Minute')
                xlim([RecordingMetrics.minBinsSPM(1) RecordingMetrics.minBinsSPM(end)])
                set(gca,'FontSize',plotProps.FigureFontSize)

                AnalysisRegionEntries = zeros(1,2); regionCounter = 1;
                RegionLabels = {}; %initialize list of region descriptions
                while doneDividing ~=1 %while there are still data regions to select...

                    disp('NOTE:  The smaller the region sizes, the less time the algorithm will require to run.  Long regions can take several hours to process.')
                    StartMin_i = input('     Please enter the x value (time) of the START of the region you want to analyze.  ');
                    EndMin_i = input('     Please enter the x value (time) of the END of the region you want to analyze.  ');


                    %check that the current region does not overlap with any previous region
                    if regionCounter>1 %if there are previously selected regions...
                        overlapTest = zeros(regionCounter-1,2); %initialize overlap test array (col 1 indicates if current region starts in a previous region, col 2 indicates if current region ends in a previous region)
                        for olTest = 1:regionCounter-1%for each previous region
                            startOL = (AnalysisRegionEntries(olTest,1) <=  StartMin_i & StartMin_i <= AnalysisRegionEntries(olTest,2)); %see if the current region starts in a previous region
                            endOL = (AnalysisRegionEntries(olTest,1) <=  EndMin_i & EndMin_i <= AnalysisRegionEntries(olTest,2)); %see if the current region ends in a previous region
                            overlapTest(olTest,:) = [startOL endOL]; %save the logicals in the overlap array
                        end %end check for current vs previous region overlap

                    else %otherwise, this is the first region being analyzed
                        overlapTest = 0; %since this is the first region, no overlap with previous regions is possible because previous regions do not exist
                    end

                    %check result of overlap test
                    if any(overlapTest(:)) %if current region overlaps with a previous region
                        disp(' ')
                        disp('Restarting selection of this region because the currently selected region overlaps with a previously selected region.  Please ensure analysis regions do not overlap. ')
                        disp(' ')
                    else %otherwise, the new region does not overlap with previous regions, and can be used in data analysis
                        hold on
                        regStartIdx = max(1,find(RecordingMetrics.minBinsSPM<max(2,StartMin_i),1,'last'));
                        regEndIdx =  max(1,find(RecordingMetrics.minBinsSPM<max(2,EndMin_i),1,'last'));
                        h2=plot(RecordingMetrics.minBinsSPM(regStartIdx:regEndIdx),RecordingMetrics.meanSPM(regStartIdx:regEndIdx),'m','LineWidth',2);
                        hold off
                        legend([h1 h2],'Entire Recording','Selected Region(s) to Analyze')
                        xlim([RecordingMetrics.minBinsSPM(1) RecordingMetrics.minBinsSPM(end)])

                        RegionLabel = input('     Please type a label for the region (i.e. Before Treatment, 10 min After Treatment, etc...). ','s');

                        regionMistake = input('     Was there a mistake in any of the entries for this region?  Type 1 if yes, 0 if no.  If yes, do not fret.  You will re-enter the correct information for this region. ');
                        if regionMistake == 0 %if there are no mistakes, save the injury/treatment
                            RegionLabels = [RegionLabels, RegionLabel]; %save the region label in the array
                            AnalysisRegionEntries(regionCounter,:) = [StartMin_i EndMin_i]; %save the region start and end in the array
                            regionCounter = regionCounter+1; %update the number of regions

                            doneDividing= input('Are you finished selecting regions?  Enter 1 if yes, 0 if no. ');

                            %If there was a mistake, the code will just repeat the loop again without saving the incorrect information
                        end


                    end %end overlap check

                end %end loop through region selection
                close(gcf) %close the region selection figure

                %Ensure that region order is chronological
                [AnalysisRegionEntries, sortInds] = sortrows(AnalysisRegionEntries,1); %sort the region entries so first entry = earliest in the recording
                RegionLabels = RegionLabels(sortInds); %sort the region labels too

                AnalysisRegions.Minutes = AnalysisRegionEntries;
                AnalysisRegions.Seconds = AnalysisRegionEntries.*60;
                AnalysisRegions.RegionLabels = RegionLabels;
                AnalysisRegions.numRegions = size(AnalysisRegions.Minutes,1); %get the number of regions to be analyzed

                save([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'AnalysisRegions') %save the regions

            else %otherwise, just seperate the recording into discrete regions
                %Note:  large regions are more likely to crash the computer because such regions might contain too many events
                %The bin size (x) is set by the user in the main function via the parameter plotProps.automatedRegionBinLength
                %The option of no manual selection bins the recording into x min sections

                if rem( max(RecordingMetrics.minBinsSPM),plotProps.automatedRegionBinLength) ~= 0 %if the end of the recording is not a full bin
                    endBinRound = max(RecordingMetrics.minBinsSPM);
                else %otherwise, the recording is divisable by the number of automatically assigned regions (without any remainder)
                    endBinRound = (plotProps.automatedRegionBinLength-rem( max(RecordingMetrics.minBinsSPM),plotProps.automatedRegionBinLength))+ max(RecordingMetrics.minBinsSPM);
                end
                binRanges = 0:plotProps.automatedRegionBinLength:endBinRound;
                binStartz = binRanges(1:end-1); 
                binEndz = binRanges(2:end)-(plotProps.automatedRegionBinLength/1000); %shift bin ends a bit so regions do not overlap
                binEndz(end) = min(binEndz(end),ProcessedData.ExperimentEndTime/60);

                %Make sure the last region is not too short, and roughly matches the length of the others (exclude if too short because it will otherwise be difficult to compare against previous regions)
                if (binEndz(end)-binStartz(end))/plotProps.automatedRegionBinLength < 0.8 %if the last region is < 80% of the bin size, exclude it from analysis
                    binStartz = binStartz(1:end-1); %exclude the very end of the recording from analysis because it is too short
                    binEndz = binEndz(1:end-1); %exclude the very end of the recording from analysis because it is too short
                end

                AnalysisRegions.Minutes = [binStartz', binEndz'];

                %find whether there are any regions with zero events
                idx = cellfun(@(y) RecordingMetrics.minBinsSPM>=y(1) & RecordingMetrics.minBinsSPM<=y(2),num2cell(AnalysisRegions.Minutes,2),'UniformOutput',false); %identify in which minute raster events fall
                emptyRegions = cell2mat(cellfun(@(x) sum(RecordingMetrics.meanSPM(x)>0)/length(RecordingMetrics.meanSPM(x)), idx,'UniformOutput',false)); 
                emptyRegions = emptyRegions<1; %create logical array if there are 0 raster events in the region

                AnalysisRegions.Minutes(emptyRegions,:) = []; %remove regions with 0 events

                AnalysisRegions.Seconds = AnalysisRegions.Minutes.*60;

                AnalysisRegions.numRegions = size(AnalysisRegions.Minutes,1); %get the number of regions to be analyzed

                save([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'AnalysisRegions') %save the array

            end %end check whether to select regions by hand or just divide the recording into regularly sized bins

        elseif AddToRegionsIndicator == 1 %otherwise, user wants to add to previously-selected regions

            %Load previously-selected regions
            load([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'AnalysisRegions');

            %Display the regions that are already selected
            figure()
            hold on
            h1=plot(RecordingMetrics.minBinsSPM,RecordingMetrics.meanSPM,'k');
            for rr = 1:AnalysisRegions.numRegions %cycle through the regions already there...
                AnalysisRegions.Minutes(rr,2)
                regStartIdx = max(1,find(RecordingMetrics.minBinsSPM<max(1,AnalysisRegions.Minutes(rr,1)),1,'last'));
                regEndIdx =  max(1,find(RecordingMetrics.minBinsSPM<max(2,AnalysisRegions.Minutes(rr,2)),1,'last'));
                h2=plot(RecordingMetrics.minBinsSPM(regStartIdx:regEndIdx),RecordingMetrics.meanSPM(regStartIdx:regEndIdx),'m','LineWidth',2);
            end
            hold off
            legend([h1 h2],'Entire Recording','Selected Region(s) to Analyze')
            xlabel('Time (min)')
            ylabel('Events per Minute After Injury or Treatment')
            xlabel('Time (min)')
            ylabel('Events per Minute')
            xlim([RecordingMetrics.minBinsSPM(1) RecordingMetrics.minBinsSPM(end)])
            set(gca,'FontSize',plotProps.FigureFontSize)

            %Add to the selection
            doneDividing = 0; %make a switch that will set to 1 once you have finished selecting regions of the data to analyze

            AnalysisRegionEntries = zeros(1,2); regionCounter = 1;
            RegionLabels = {}; %initialize list of region descriptions
            while doneDividing ~=1 %while there are still data regions to select...

                disp('NOTE:  The smaller the region sizes, the less time the algorithm will require to run.  Long regions can take several hours to process.')
                StartMin_i = input('     Please enter the x value (time) of the START of the region you want to analyze.  ');
                EndMin_i = input('     Please enter the x value (time) of the END of the region you want to analyze.  ');


                %check that the current region does not overlap with any previous region
                if regionCounter>1 %if there are previously selected regions...
                    overlapTest = zeros(regionCounter-1,2); %initialize overlap test array (col 1 indicates if current region starts in a previous region, col 2 indicates if current region ends in a previous region)
                    for olTest = 1:regionCounter-1%for each previous region
                        startOL = (AnalysisRegionEntries(olTest,1) <=  StartMin_i & StartMin_i <= AnalysisRegionEntries(olTest,2)); %see if the current region starts in a previous region
                        endOL = (AnalysisRegionEntries(olTest,1) <=  EndMin_i & EndMin_i <= AnalysisRegionEntries(olTest,2)); %see if the current region ends in a previous region
                        overlapTest(olTest,:) = [startOL endOL]; %save the logicals in the overlap array
                    end %end check for current vs previous region overlap

                else %otherwise, this is the first region being analyzed
                    overlapTest = 0; %since this is the first region, no overlap with previous regions is possible because previous regions do not exist
                end

                %check result of overlap test
                if any(overlapTest(:)) %if current region overlaps with a previous region
                    disp(' ')
                    disp('Restarting selection of this region because the currently selected region overlaps with a previously selected region.  Please ensure analysis regions do not overlap. ')
                    disp(' ')
                else %otherwise, the new region does not overlap with previous regions, and can be used in data analysis
                    hold on
                    regStartIdx = max(1,find(RecordingMetrics.minBinsSPM<max(2,StartMin_i),1,'last'));
                    regEndIdx =  max(1,find(RecordingMetrics.minBinsSPM<max(2,EndMin_i),1,'last'));
                    h2=plot(RecordingMetrics.minBinsSPM(regStartIdx:regEndIdx),RecordingMetrics.meanSPM(regStartIdx:regEndIdx),'m','LineWidth',2);
                    hold off
                    legend([h1 h2],'Entire Recording','Selected Region(s) to Analyze')
                    xlim([RecordingMetrics.minBinsSPM(1) RecordingMetrics.minBinsSPM(end)])

                    RegionLabel = input('     Please type a label for the region (i.e. Before Treatment, 10 min After Treatment, etc...). ','s');

                    regionMistake = input('     Was there a mistake in any of the entries for this region?  Type 1 if yes, 0 if no.  If yes, do not fret.  You will re-enter the correct information for this region. ');
                    if regionMistake == 0 %if there are no mistakes, save the injury/treatment
                        RegionLabels = [RegionLabels, RegionLabel]; %save the region label in the array
                        AnalysisRegionEntries(regionCounter,:) = [StartMin_i EndMin_i]; %save the region start and end in the array
                        regionCounter = regionCounter+1; %update the number of regions

                        doneDividing= input('Are you finished selecting regions?  Enter 1 if yes, 0 if no. ');

                        %If there was a mistake, the code will just repeat the loop again without saving the incorrect information
                    end


                end %end check whether to add to the analysis regions or whether to select them fresh

            end %end check whether still selecting regions

            %combine new and old regions
            AnalysisRegionEntries = [AnalysisRegions.Minutes; AnalysisRegionEntries];
            RegionLabels = [AnalysisRegions.RegionLabels, RegionLabels];

            %Ensure that region order is chronological
            [AnalysisRegionEntries, sortInds] = sortrows(AnalysisRegionEntries,1); %sort the region entries so first entry = earliest in the recording
            RegionLabels = RegionLabels(sortInds); %sort the region labels too

            %Now add to the saved regions and save the updated and sorted region list
            AnalysisRegions.Minutes = AnalysisRegionEntries;
            AnalysisRegions.Seconds = 60.*AnalysisRegions.Minutes;
            AnalysisRegions.RegionLabels = RegionLabels;
            AnalysisRegions.numRegions = size(AnalysisRegions.Minutes,1); %get the number of regions to be analyzed

            save([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'AnalysisRegions')

        end  %end check whether to slect new regions or load old regions

    else %otherwise just load previously-set analysis regions
        load([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'AnalysisRegions');

    end %end check whether to set or load analysis regions

else %otherwise, just load previously-selected regions and injuries/treatments

    load([pwd '\Output Files\' plxFileName '\InjuryIndicies.mat'],'InjuryIndicies');
    load([pwd '\Output Files\' plxFileName '\AnalysisRegions.mat'],'AnalysisRegions');

end %end check whether to label injuries/treatments and analysis regions, or load previously-selected


%% Now plot the results
%Plot the mean event count/minute with injuries/treatments labeled

spikesPerMin = RecordingMetrics.individualSPM; %get individual signal events/min for plotting
meanSPM = RecordingMetrics.meanSPM; %get mean events/min for plotting
stdSPM = RecordingMetrics.stdSPM; %get the standard deviation in events/min
minBins = RecordingMetrics.minBins; %get each min

ymax = 1.1*max([max(cellfun(@max,spikesPerMin)) , max(meanSPM+stdSPM)]);
ymin = min([min(cellfun(@min,spikesPerMin)) , min(meanSPM-stdSPM), -1]);


figure()
hold on
%Plot injury or treatment if applicable
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury/treatment
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        if InjuryIndicies.InjuryStartMin(j) > 0 %if the injury/treatment occurred during the recording
            injuryStartIndOnPlot = find(RecordingMetrics.minBinsSPM<InjuryIndicies.InjuryStartMin(j),1,'last')+1;
            injuryEndIndOnPlot = find(RecordingMetrics.minBinsSPM>InjuryIndicies.InjuryEndMin(j),1,'first')-1;
            if injuryStartIndOnPlot ~= injuryEndIndOnPlot
                x=[injuryStartIndOnPlot injuryEndIndOnPlot injuryEndIndOnPlot injuryStartIndOnPlot];
                A=fill(x,[ymin-10 ymin-10 ymax*1.2 ymax*1.2],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none');
            else
                x=[injuryStartIndOnPlot-0.1 injuryEndIndOnPlot+0.1 injuryEndIndOnPlot+0.1 injuryStartIndOnPlot-0.1];
                A=fill([injuryStartIndOnPlot-0.1 injuryEndIndOnPlot+0.1 injuryEndIndOnPlot+0.1 injuryStartIndOnPlot-0.1],[ymin-10 ymin-10 ymax*1.2 ymax*1.2],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none');

            end
            hh=hatchfill(A,'cross',44,10,[0.6 0.6 0.6]);
            set(hh,'color',[0.7 0.7 0.7]);
            plot([injuryStartIndOnPlot injuryStartIndOnPlot],[ymin-10 ymax*1.2],'--k','LineWidth',plotProps.lineWidth/2)
            plot([injuryEndIndOnPlot injuryEndIndOnPlot],[ymin-10 ymax*1.2],'--k','LineWidth',plotProps.lineWidth/2)

            text(injuryStartIndOnPlot+(injuryEndIndOnPlot-injuryStartIndOnPlot)/2,ymax*.9,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            
        end
    end %end cycle through injuries/treatments

end %end check whether to plot an injury or treatment

p = cellfun(@plot,spikesPerMin); %plot rasters (y axis is cell, x axis is time)
for i = 1:length(RecordingMetrics.individualSPM) %cycle through cells...
    p(i).Marker = '.';
    p(i).LineStyle = '-';
    p(i).Color = [0.7 0.7 1];
    p(i).LineWidth = plotProps.lineWidth/2;
    p(i).MarkerSize = plotProps.markerSize/2;
end

x=[1:length(meanSPM) fliplr(1:length(meanSPM))];
y=[meanSPM+stdSPM fliplr(meanSPM-stdSPM)];
fill(x,y,'b','FaceAlpha',plotProps.faceAlpha)
plot(meanSPM,'.-b','LineWidth',plotProps.lineWidth,'MarkerSize',plotProps.markerSize)
plot(meanSPM-stdSPM,'--b','LineWidth',plotProps.lineWidth,'MarkerSize',plotProps.markerSize)
plot(meanSPM+stdSPM,'--b','LineWidth',plotProps.lineWidth,'MarkerSize',plotProps.markerSize)
h1 = plot(NaN,NaN,'.-','Color',[0.7 0.7 1],'LineWidth',plotProps.lineWidth/2,'MarkerSize',plotProps.markerSize/2);
h2 = plot(NaN,NaN,'.-b','LineWidth',plotProps.lineWidth,'MarkerSize',plotProps.markerSize);
hold off
legend([h1,h2],'Individual Signals','Mean')
set(gca,'XTick',minBins(2): (minBins(end)-minBins(2))/10: minBins(end),'XTickLabel',minBins(2): (minBins(end)-minBins(2))/10: minBins(end))
xlabel('Time (min)')
ylabel('Events per Minute')
xlim([minBins(2) minBins(end)])
ylim([ymin-1 ymax+1])
box on
title('Changes in Events/Minute Through Time')
set(gca,'FontSize',plotProps.FigureFontSize)
grid on



%Plot the regions being analyzed
numRegions =  size(AnalysisRegions.Minutes,1);

if numRegions > size(plotProps.regionColors,1) %if your colormap doesn't have enough different colors to plot a unique color for all regions...
    plotProps.regionColors = colorcube(numRegions+1);
end

figure()
hold on
%Plot injury or treatment if applicable
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury/treatment
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        if InjuryIndicies.InjuryStartMin(j) > 0 %if the injury/treatment occurred during the recording
            injuryStartIndOnPlot = find(RecordingMetrics.minBinsSPM<InjuryIndicies.InjuryStartMin(j),1,'last')+1;
            injuryEndIndOnPlot = find(RecordingMetrics.minBinsSPM>InjuryIndicies.InjuryEndMin(j),1,'first')-1;
            if injuryStartIndOnPlot ~= injuryEndIndOnPlot
                x=[injuryStartIndOnPlot injuryEndIndOnPlot injuryEndIndOnPlot injuryStartIndOnPlot];
                A=fill(x,[ymin-10 ymin-10 ymax*1.2 ymax*1.2],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none');
            else
                x=[injuryStartIndOnPlot-0.1 injuryEndIndOnPlot+0.1 injuryEndIndOnPlot+0.1 injuryStartIndOnPlot-0.1];
                A=fill([injuryStartIndOnPlot-0.1 injuryEndIndOnPlot+0.1 injuryEndIndOnPlot+0.1 injuryStartIndOnPlot-0.1],[ymin-10 ymin-10 ymax*1.2 ymax*1.2],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none');

            end
            hh=hatchfill(A,'cross',44,10,[0.6 0.6 0.6]);
            set(hh,'color',[0.7 0.7 0.7]);
            plot([injuryStartIndOnPlot injuryStartIndOnPlot],[ymin-10 ymax*1.2],'--k','LineWidth',plotProps.lineWidth/2)
            plot([injuryEndIndOnPlot injuryEndIndOnPlot],[ymin-10 ymax*1.2],'--k','LineWidth',plotProps.lineWidth/2)
        end
    end %end cycle through injuries/treatments

end %end check whether to plot an injury or treatment

plot(RecordingMetrics.minBinsSPM,RecordingMetrics.meanSPM,'k','LineWidth',plotProps.lineWidth)

%Now plot analysis regions...

%first get the font size to label analysis regions
if AnalysisRegions.numRegions < 100
    textSz = plotProps.FigureFontSize*0.8;
else %if there are many regions to label, the text box font size needs to decrease accordingly
    textSz = plotProps.FigureFontSize*0.4;
end

for i = 1:AnalysisRegions.numRegions
    Inds = max(1,AnalysisRegions.Minutes(i,:)); %maxe sure index >= 1
    Inds = min(length(RecordingMetrics.minBinsSPM),Inds); %make sure index <= total recording length
    Inds = round(Inds); %make sure values are integers for plotting
    fill([Inds(1) Inds(2) Inds(2) Inds(1)],[0 0 max(RecordingMetrics.meanSPM).*1.4 max(RecordingMetrics.meanSPM).*1.4],plotProps.regionColors(i,:),'facealpha',0.2,'LineStyle','none');
    plot([Inds(1) Inds(1)], [0 max(RecordingMetrics.meanSPM).*1.4]  ,'--','Color',plotProps.regionColors(i,:).*0.85,'LineWidth',plotProps.lineWidth/2)
    plot([Inds(2) Inds(2)], [0 max(RecordingMetrics.meanSPM).*1.4]  ,'--','Color',plotProps.regionColors(i,:).*0.85,'LineWidth',plotProps.lineWidth/2)
    plot(RecordingMetrics.minBinsSPM(Inds(1):Inds(2)),RecordingMetrics.meanSPM(Inds(1):Inds(2)),'Color',plotProps.regionColors(i,:).*0.7,'LineWidth',plotProps.lineWidth)
    if ~isfield(AnalysisRegions,'RegionLabels') %if user automated region selection and hence did not type region labels...
        if AnalysisRegions.numRegions < 30 %if you have few enough regions to write the full label
            text(Inds(1)+(Inds(2)-Inds(1))/2,ceil(max(RecordingMetrics.meanSPM(Inds(1):Inds(2))).*1.2),['Region ' num2str(i)],'FontSize',textSz,'HorizontalAlignment', 'center','Color',plotProps.regionColors(i,:).*0.5,'FontAngle','italic')
        else %otherwise, so many regions that writing the full labels will be illegible, so just write the region number
            text(Inds(1)+(Inds(2)-Inds(1))/2,ceil(max(RecordingMetrics.meanSPM(Inds(1):Inds(2))).*1.2), num2str(i),'FontSize',textSz,'HorizontalAlignment', 'center','Color',plotProps.regionColors(i,:).*0.5,'FontAngle','italic')
        end
    else %otherwise user assigned specific region labels
        text(Inds(1)+(Inds(2)-Inds(1))/2,ceil(max(RecordingMetrics.meanSPM(Inds(1):Inds(2))).*1.2), AnalysisRegions.RegionLabels{i},'FontSize',textSz,'HorizontalAlignment', 'center','Color',plotProps.regionColors(i,:).*0.5,'FontAngle','italic')
    end %end check whether user assigned specific region labels
end

%Label treatments (the second for loop here is so the treatment text shows up on top of the plots.  If it is above, the text boxes can get buried under the other parts of the plot)
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury/treatment
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        if InjuryIndicies.InjuryStartMin(j) > 0 %if the injury/treatment occurred during the recording
            injuryStartIndOnPlot = find(RecordingMetrics.minBinsSPM<InjuryIndicies.InjuryStartMin(j),1,'last')+1;
            injuryEndIndOnPlot = find(RecordingMetrics.minBinsSPM>InjuryIndicies.InjuryEndMin(j),1,'first')-1;
            
            text(injuryStartIndOnPlot+(injuryEndIndOnPlot-injuryStartIndOnPlot)/2,max(RecordingMetrics.meanSPM).*1.3,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')

        end
    end %end cycle through injuries/treatments

end %end check whether to plot an injury or treatment

hold off
xlabel('Time (min)')
ylabel('Average Events per Minute')
title('Regions Selected for Analysis')
xlim([0.99*RecordingMetrics.minBinsSPM(1) RecordingMetrics.minBinsSPM(end)+1])
ylim([0 max(RecordingMetrics.meanSPM).*1.4])
box on
set(gca,'FontSize',plotProps.FigureFontSize)
grid on

end %end function