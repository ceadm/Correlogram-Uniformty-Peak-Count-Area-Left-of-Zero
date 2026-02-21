function ProcessedData = plotStatsInEachRegion(AnalysisRegions,InjuryIndicies,ProcessedData,plotProps,sw) 
%This function plots the firing count and interval histograms either by 
%analysis region or before vs. after treatment based on user preference
%(instead of for the entire recording)

%% Get histogram bins
%Event count bins for histogram (matches the overall recording bins)
SCbins = linspace(-4,1+log10(max(ProcessedData.FiringCountOfEachCell)),round(ProcessedData.CellCount/1.5)); SCbins = 10.^SCbins;

%Interevent interval bins for histogram (matches the overall recording bins)
SIbins = linspace(-5,4,100); SIbins = 10.^SIbins;

%Save in the data structure
ProcessedData.RegionStats.SCbins = SCbins;
ProcessedData.RegionStats.SIbins = SIbins;

%If the user says they wish to seperate the recording based on treatments, not analysis regions, but there is no treatment during the recording, plot by analysis region
if sw.plotRecordingStatsForEachDataRegion == 0 && isempty(find(InjuryIndicies.InjuryStartSec>0,1,'first')); sw.plotRecordingStatsForEachDataRegion = 1; end

%Now check how to plot
if sw.plotRecordingStatsForEachDataRegion == 1 %if plotting stats by analysis region (independent of treatments)

    %Get the number of rows and columns for the subplots of correlograms in each data analysis region if there are so many regions that tiled figure layouts will be used
    c = min(10,AnalysisRegions.numRegions); %number of columns for tile()
    r = ceil(AnalysisRegions.numRegions/c); %number of rows for tile()

    %% Label (or fetch labels) for each recording region
    if ~isfield(AnalysisRegions,'RegionLabels') %if user automated region selection and hence did not type region labels...
        groupOrder = cell(1,AnalysisRegions.numRegions); %get each analysis region
        for rr = 1:AnalysisRegions.numRegions %for each region...
            groupOrder{rr} = 'Region' + " " + num2str(rr); %get the region name
        end %end cycle through regions

    else %otherwise user assigned specific region labels, so use those instead
        groupOrder = AnalysisRegions.RegionLabels;
    end

    %Save the groupings and times within each group in the ProcessedData structure
    ProcessedData.RegionStats.groupOrder = groupOrder;
    ProcessedData.RegionStats.timeRange = AnalysisRegions.Seconds';
    timeRange = AnalysisRegions.Seconds';

    %% Make Figures

    %Initialize figures
    if AnalysisRegions.numRegions <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
        f_spikeCount = figure;
        f_interspikeInterval = figure;
    else %otherwise, there are so many analysis regions that they shouldn't all be on the same plot
        f_spikeCount = figure;
        sct = tiledlayout(r,c,'Padding','none','TileSpacing','compact');
        f_interspikeInterval = figure;
        isit = tiledlayout(r,c,'Padding','none','TileSpacing','compact');
    end

    %Now plot the histograms for each region
    for rr = 1:AnalysisRegions.numRegions %for each analysis region...

        reg_field = strcat('Region',num2str(rr)); %get the region name to save in the AnalysisRegions data structure

        %Count the number of events that occurred within the given analysis region
        spikesInRegion = cellfun(@(x) x>=AnalysisRegions.Seconds(rr,1) & x<=AnalysisRegions.Seconds(rr,2), ProcessedData.Raster,'UniformOutput',false); %Find spikes that occurred within the given analysis region
        regionSpikeCount = cell2mat(cellfun(@(x) sum(x), spikesInRegion,'UniformOutput',false));  %Now count the spikes that occured within the given analysis region
        regionSpikeCount = regionSpikeCount/(timeRange(2,rr)-timeRange(1,rr)); %convert to spikes per minute to account for regions of varying length
        fractionOfSpikesThatOccurrInThisRegion = regionSpikeCount./ProcessedData.FiringCountOfEachCell; %Get the fraction of the cell's total event count that occurred in the given analysis region

        %Indicies to look at for interevent intervals
        INDS_interspikeIntervalsInRegion = cellfun(@(x) max(1,find(x)-1), spikesInRegion,'UniformOutput',false); %event interval indicies
        interspikeIntervalsInRegion = cellfun(@(x,y) x(y),  ProcessedData.SpikeIntervals,INDS_interspikeIntervalsInRegion,'UniformOutput',false); %list of event intervals in the given region
        
        
        %Save the stats in the data structure
        ProcessedData.RegionStats.(reg_field).spikesInRegion = spikesInRegion;
        ProcessedData.RegionStats.(reg_field).regionSpikeCount = regionSpikeCount;
        ProcessedData.RegionStats.(reg_field).fractionOfSpikesThatOccurrInThisRegion =  fractionOfSpikesThatOccurrInThisRegion;
        ProcessedData.RegionStats.(reg_field).INDS_interspikeIntervalsInRegion = INDS_interspikeIntervalsInRegion;
        ProcessedData.RegionStats.(reg_field).interspikeIntervalsInRegion = interspikeIntervalsInRegion;

        interspikeIntervalsInRegion = cell2mat(interspikeIntervalsInRegion'); %convert interval cell array to a numerical vector for plotting

        %Plot event count of each cell
        if AnalysisRegions.numRegions <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
            figure(f_spikeCount)
            hold on
            histogram(regionSpikeCount,SCbins,'EdgeColor',plotProps.regionColors(rr,:),'FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
            hold off
        else %otherwise, there are so many analysis regions that they shouldn't all be on the same plot
            ax=nexttile(sct);
            histogram(ax,regionSpikeCount,SCbins,'EdgeColor',plotProps.regionColors(rr,:),'FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
            title(groupOrder{rr})
        end
        box on
        grid on
        set(gca,'XScale','log')
        set(gca,'FontSize',plotProps.FigureFontSize)
        set(gca,'XTick',10.^[-4:1:round(log10(max(SCbins)))])
        xlim(10.^[-4 round(log10(max(SCbins)))])

        %Now plot interevent interval
        if AnalysisRegions.numRegions <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
            figure(f_interspikeInterval)
            hold on
            histogram(interspikeIntervalsInRegion,SIbins,'Normalization','probability','EdgeColor',plotProps.regionColors(rr,:),'FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
            hold off
        else %otherwise, there are so many analysis regions that they shouldn't all be on the same plot
            ax=nexttile(isit);
            histogram(ax,interspikeIntervalsInRegion,SIbins,'Normalization','probability','EdgeColor',plotProps.regionColors(rr,:),'FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
            title(groupOrder{rr},'FontSize',plotProps.FigureFontSize)
        end
        %set(gca,'YScale','log')
        set(gca,'XScale','log')
        box on
        grid on
        set(gca,'FontSize',plotProps.FigureFontSize)
        set(gca,'XTick',10.^[round(log10(min(SIbins))-1):1:round(log10(max(SIbins)))])
        xlim(10.^[round(log10(min(SIbins))-1) round(log10(max(SIbins)))])

    end %end cycle through regions

    %Label plot axes
    if AnalysisRegions.numRegions <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
        figure(f_spikeCount)
        legend(groupOrder)
        xlabel('Number of Times a Signal Fired (Count / Min)','FontSize',plotProps.FigureFontSize)
        ylabel('Number of Signals','FontSize',plotProps.FigureFontSize)
        title('Firing Count Histogram','FontSize',plotProps.FigureFontSize)

        figure(f_interspikeInterval)
        legend(groupOrder)
        xlabel('Interval Between Events (s)')
        ylabel('Probability')
        title('Event Interval Histogram')

    else %otherwise, there are so many analysis regions that they shouldn't all be on the same plot
        xlabel(sct,'Number of Times a Signal Fired (Count / Min)','FontSize',plotProps.FigureFontSize);
        ylabel(sct,'Number of Signals','FontSize',plotProps.FigureFontSize);

        xlabel(isit,'Interval Between Events (s)','FontSize',plotProps.FigureFontSize);
        ylabel(isit,'Probability','FontSize',plotProps.FigureFontSize);
    end %end check whether all regions are on same plot or not



else %otherwise, divide the data according to treatments

    if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there are injuries or treatments...

        %if treatments were administered before the recording started, do not include them to split up the recording
        b =  find(InjuryIndicies.InjuryStartSec>0,1,'first'); %find treatments administered during the recording

        if ~isempty(b) %if there were any treatments during recording, divide up the recording

            subtractInd = b-1; %Figure how many treatments were applied before the recording

            groupOrder = cell(1,1+InjuryIndicies.NumberOfInjuriesOrTreatments-(b-1)); %initialize array to name each region of the recording
            timeRange = zeros(2,1+InjuryIndicies.NumberOfInjuriesOrTreatments-(b-1)); %initialize array to store times of the region in the recording, row 1 = each region start, row 2 = each region end
            for rr = b:1+InjuryIndicies.NumberOfInjuriesOrTreatments %for each region...

                %Note any treatments administered before the recording with an extra string
                if rr == b && b>1 %if region is between the last treatment before the recording and the first treatment during recording, say so with the string
                    preStartStr = ' (administered before recording) ';
                else %otherwise not comparing with something administered before recording
                    preStartStr = ' ';
                end

                %Name the region and get the region start and end (in seconds)
                if rr == 1
                    groupOrder{rr-subtractInd} = 'Before' + " " + InjuryIndicies.InjuryLabels{rr};
                    timeRange(1,rr-subtractInd) = 0; %start time is the recording start (t=0)
                    timeRange(2,rr-subtractInd) = InjuryIndicies.InjuryStartSec(rr)-0.001; %end time is 1 ms before the start of the first treatment
                elseif rr>1 && rr<=InjuryIndicies.NumberOfInjuriesOrTreatments
                    groupOrder{rr-subtractInd} = 'Between' + " " + InjuryIndicies.InjuryLabels{rr-1} + preStartStr + "and " + InjuryIndicies.InjuryLabels{rr};
                    timeRange(1,rr-subtractInd) = max(0,0.001+InjuryIndicies.InjuryEndSec(rr-1)); %start time is 1 ms after the end of the previous treatment or 0 (if previous treatment occurred before recording start
                    timeRange(2,rr-subtractInd) = InjuryIndicies.InjuryStartSec(rr)-0.001; %end time is 1 ms before the start of the next treatment
                else %otherwise, last region
                    groupOrder{rr-subtractInd} = 'After' + " " + InjuryIndicies.InjuryLabels{rr-1};
                    timeRange(1,rr-subtractInd) = max(0,0.001+InjuryIndicies.InjuryEndSec(rr-1)); %start time is 1 ms after the end of the previous treatment or 0 (if previous treatment occurred before recording start
                    timeRange(2,rr-subtractInd) = ProcessedData.ExperimentEndTime; %end time is the recording end
                end

            end %end cycle through regions

            %Save the groupings and times within each group in the ProcessedData structure
            ProcessedData.RegionStats.groupOrder = groupOrder;
            ProcessedData.RegionStats.timeRange = timeRange;

            %If there are so many treatments that tiled figures are needed, get tiled dimensions
            if size(groupOrder,2) > plotProps.numberOfRegionsForASinglePlot
                %Get the number of rows and columns for the subplots of histograms in each data analysis region if there are so many regions that tiled figure layouts will be used
                c = min(10,size(groupOrder,2)); %number of columns for tile()
                r = ceil(size(groupOrder,2)/c); %number of rows for tile()
            end %end check whether plotting histograms all on the same plot or in a tiled layout

            %Initialize figures
            if size(groupOrder,2) <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
                f_spikeCount = figure;
                f_interspikeInterval = figure;
            else %otherwise, there are so many analysis regions that they shouldn't all be on the same plot
                f_spikeCount = figure;
                sct = tiledlayout(r,c,'Padding','none','TileSpacing','compact');
                f_interspikeInterval = figure;
                isit = tiledlayout(r,c,'Padding','none','TileSpacing','compact');
            end

            %Now plot the histograms for each region
            for rr = 1:size(groupOrder,2) %for each analysis region...

                reg_field = strcat('Region',num2str(rr)); %get the region name to save in the AnalysisRegions data structure

                %Count the number of spikes that occurred within the given analysis region
                spikesInRegion = cellfun(@(x) x>=timeRange(1,rr) & x<=timeRange(2,rr), ProcessedData.Raster,'UniformOutput',false); %Find spikes that occurred within the given analysis region
                regionSpikeCount = cell2mat(cellfun(@(x) sum(x), spikesInRegion,'UniformOutput',false));  %Now count the spikes that occured within the given analysis region
                regionSpikeCount = regionSpikeCount/(timeRange(2,rr)-timeRange(1,rr)); %convert to spikes per minute to account for regions of varying length
                fractionOfSpikesThatOccurrInThisRegion = regionSpikeCount./ProcessedData.FiringCountOfEachCell; %Get the fraction of the cell's total event count that occurred in the given analysis region

                %Indicies to look at for interevent intervals
                INDS_interspikeIntervalsInRegion = cellfun(@(x) max(1,find(x)-1), spikesInRegion,'UniformOutput',false); %event interval indicies
                interspikeIntervalsInRegion = cellfun(@(x,y) x(y),  ProcessedData.SpikeIntervals,INDS_interspikeIntervalsInRegion,'UniformOutput',false); %list of event intervals in the given region

                %Save the stats in the data structure
                ProcessedData.RegionStats.(reg_field).spikesInRegion = spikesInRegion;
                ProcessedData.RegionStats.(reg_field).regionSpikeCount = regionSpikeCount;
                ProcessedData.RegionStats.(reg_field).fractionOfSpikesThatOccurrInThisRegion =  fractionOfSpikesThatOccurrInThisRegion;
                ProcessedData.RegionStats.(reg_field).INDS_interspikeIntervalsInRegion = INDS_interspikeIntervalsInRegion;
                ProcessedData.RegionStats.(reg_field).interspikeIntervalsInRegion = interspikeIntervalsInRegion;


                interspikeIntervalsInRegion = cell2mat(interspikeIntervalsInRegion'); %convert interval cell array to a numerical vector for plotting

                %Plot event count of each cell
                if size(groupOrder,2) <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
                    regionColors = jet(size(groupOrder,2));
                    figure(f_spikeCount)
                    hold on
                    histogram(regionSpikeCount,SCbins,'EdgeColor',regionColors(rr,:),'FaceColor',regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
                    hold off
                else %otherwise, there are so many analysis regions that they shouldn't all be on the same plot
                    regionColors = jet(size(groupOrder,2));
                    ax=nexttile(sct);
                    histogram(ax,regionSpikeCount,SCbins,'EdgeColor',regionColors(rr,:),'FaceColor',regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
                    title(groupOrder{rr})
                end
                box on
                grid on
                set(gca,'XScale','log')
                set(gca,'FontSize',plotProps.FigureFontSize)
                set(gca,'XTick',10.^[-4:1:round(log10(max(SCbins)))])
                xlim(10.^[-4 round(log10(max(SCbins)))])

                %Now plot interevent interval
                if size(groupOrder,2) <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
                    figure(f_interspikeInterval)
                    hold on
                    histogram(interspikeIntervalsInRegion,SIbins,'Normalization','probability','EdgeColor',regionColors(rr,:),'FaceColor',regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
                    hold off
                else %otherwise, there are so many snalysis regions that they shouldn't all be on the same plot
                    ax=nexttile(isit);
                    histogram(ax,interspikeIntervalsInRegion,SIbins,'Normalization','probability','EdgeColor',regionColors(rr,:),'FaceColor',regionColors(rr,:),'FaceAlpha',max(0.3,plotProps.faceAlpha))
                    title(groupOrder{rr},'FontSize',plotProps.FigureFontSize)
                end
                %set(gca,'YScale','log')
                set(gca,'XScale','log')
                box on
                grid on
                set(gca,'FontSize',plotProps.FigureFontSize)
                set(gca,'XTick',10.^[round(log10(min(SIbins))-1):1:round(log10(max(SIbins)))])
                xlim(10.^[round(log10(min(SIbins))-1) round(log10(max(SIbins)))])

            end %end cycle through regions

            if size(groupOrder,2) <= plotProps.numberOfRegionsForASinglePlot %if there are few enough regions to plot them all on the same plot...
                figure(f_spikeCount)
                legend(groupOrder)
                xlabel('Number of Times a Signal Fired (Count / Min)','FontSize',plotProps.FigureFontSize)
                ylabel('Number of Signals','FontSize',plotProps.FigureFontSize)
                title('Firing Count Histogram','FontSize',plotProps.FigureFontSize)

                figure(f_interspikeInterval)
                legend(groupOrder)
                xlabel('Interval Between Events (s)')
                ylabel('Probability')
                title('Event Interval Histogram')

            else %otherwise, there are so many analysis regions that they shouldn't all be on the same plot
                xlabel(sct,'Number of Times a Signal Fired (Count / Min)','FontSize',plotProps.FigureFontSize);
                ylabel(sct,'Number of Signals','FontSize',plotProps.FigureFontSize);

                xlabel(isit,'Interval Between Events (s)','FontSize',plotProps.FigureFontSize);
                ylabel(isit,'Probability','FontSize',plotProps.FigureFontSize);
            end %end check whether all regions are on same plot or not


        end %end second check whether there are injuries or treatments during the recording (or whether treatments are only before recording)

    end %end check whether there are injuries or treatments at all

end %end check whether plotting by data analysis regions or by treatment times

end %end function

