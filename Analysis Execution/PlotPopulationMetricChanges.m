function PlotPopulationMetricChanges(ProcessedData,RecordingMetrics,AnalysisRegions,InjuryIndicies,Correlograms,plotProps,sw)
%Make violin plots of correlogram metric distributions for all analyzed regions

%% Initialize arrays to make violin plots of correlogram metrics
if plotProps.includeAutocorrelogramsInStatistics == 1 %if user wants to include autocorrelograms in the violins...
    b1 = zeros(size(RecordingMetrics.CorrelogramMetrics.Region1.leaderProbMat)); %include all indicies
else %otherwise, only include crosscorrelograms in violins
    %find autocorrelograms in the matrix (same as the signal count, unless only select references are specified)
    b1=eye(size(RecordingMetrics.CorrelogramMetrics.Region1.leaderProbMat)); %find the matrix diagonal/autocorrelograms
end

%If only selecting a subset of signals as the reference signal, take only the relevent indicies (of those comparisons only)
b2=ones(ProcessedData.CellCount,ProcessedData.CellCount); b2(:,RecordingMetrics.CorrelogramMetrics.refIndicies)=0;

b=(b1 & ~b2); autoCorNum = sum(b(:));

leaderProbArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
numPeaksArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
uniformityArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
uniformitypArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);

%For peak positions (which provide frequencies for any rhythmic activity), each region can have a different number of peaks.  Therefore, regions with less than the maximum number of peaks need to be padded with NaN (here -1000 because violin() does not accept NaN) so that they are all the same size for violin plot 
MaxNumPk = 0;
for rr = 1:AnalysisRegions.numRegions %find the maximum possible number of peaks for all regions
    reg_field = strcat('Region',num2str(rr)); %convert the name of the region to a variable name to call the correct correlogram substructure
    MaxNumPk = max(MaxNumPk,max(RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks));
end
PkArray = NaN(MaxNumPk,AnalysisRegions.numRegions); 
PkTArray = NaN(MaxNumPk,AnalysisRegions.numRegions); 

groupOrderDelta = cell(1,AnalysisRegions.numRegions); %initialize name array for each region difference

if ~isfield(AnalysisRegions,'RegionLabels') %if user automated region selection and hence did not type region labels...
    groupOrder = cell(1,AnalysisRegions.numRegions); %get each analysis region
    for rr = 1:AnalysisRegions.numRegions %for each region...
        groupOrder{rr} = 'Region' + " " + num2str(rr); %get the region name
        if rr == 1 %now get the name for the difference between the present and previous region
            groupOrderDelta{rr} = 'First Region'; %Mark the first region (so no predecessor to compare against)
        else %otherwise mark the region that comes before
            groupOrderDelta{rr} = groupOrder{rr} + " " + '-' + " " + groupOrder{rr-1}; %get the region name to access the data structure;
        end
    end %end cycle through regions

else %otherwise user assigned specific region labels, so use those instead
    groupOrder = AnalysisRegions.RegionLabels;
    for rr = 1:AnalysisRegions.numRegions %for each region...
        if rr == 1 %now get the name for the difference between the present and previous region
            groupOrderDelta{rr} = 'First Region is' + " " + AnalysisRegions.RegionLabels{rr}; %Mark the first region (so no predecessor to compare against)
        else %otherwise mark the region that comes before
            groupOrderDelta{rr} = groupOrder{rr} + " " + '-' +  " " + groupOrder{rr-1}; %get the region name to access the data structure;
        end
    end %end cycle through regions
end

%% Find where the regions occur relative to injury or treatment
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury/treatment in the recording

    InjuryInds = zeros(3,InjuryIndicies.NumberOfInjuriesOrTreatments); %initialize array to save x coordinate of injuries/treatments for the violin plots (row 1) as well as the shading width (row 2) to indicate where the region occurred, and flag whether the injury occurred before the start of the recording (row 3)

    for j = 1: InjuryIndicies.NumberOfInjuriesOrTreatments %cycle through treatments/injuries

        relativeStartTimez  = AnalysisRegions.Minutes-InjuryIndicies.InjuryStartMin(j); %get the time difference between each region and the injury start
        relativeEndTimez  = AnalysisRegions.Minutes-InjuryIndicies.InjuryEndMin(j); %get the time difference between each region and the injury end

        endNegs = relativeEndTimez<0; %find where a region started (col 1) or ended (col 2) before an injury/treatment end
        endPos = relativeEndTimez>0; %find where a region started (col 1) or ended (col 2) after an injury/treatment end
        startNegs = relativeStartTimez<0; %find where a region started (col 1) or ended (col 2) before an injury/treatment end

        startNegSum = sum(startNegs,2); %if this sum is 2, then the region started and ended before the injury/treatment started.
        endNegSum = sum(endNegs,2); %if this sum is 2, then the region started and ended before the injury/treatment ended.

        endPosSum = sum(endPos,2); %if this sum is 2, then the region started and ended after the injury/treatment ended.

        simultaneousRegAndInj = endNegSum ~= startNegSum; %find whether any region was both in and out of the treatment range/simultaneous with the treatment

        if any(simultaneousRegAndInj) %if any row had both positive and negative times/occurred during an injury...
            injuryIndj = find(simultaneousRegAndInj);
            shadeWidthj = 0.4;
            startFlagj = 0; %flag=1 only if the injury occurred before the recording start

        else %otherwise the regions are either completely before or completely after the given injury/treatment
            locStartBetween = find(endNegSum==0,1,'first'); %find the first region after the injury/treatment
            locEndBetween = find(endPosSum==0,1,'last'); %find the last region before the injury/treatment

            if ~isempty(locEndBetween) && ~isempty(locStartBetween) %if the injury/treatment was during the recording
                injuryIndj = locStartBetween-(locStartBetween - locEndBetween)/2; %the index of the injury/treatment is between these two regions
                shadeWidthj = 0.03;
                startFlagj = 0;  %flag=1 if the injury occurred before the recording start

            elseif isempty(locEndBetween) %if the treatment/injury was before the recording started or after all selected regions
                if InjuryIndicies.InjuryStartMin(j) > 0 %injury/treatment was after all analyzed regions
                    injuryIndj = locStartBetween-.5; %the index of the injury/treatment is between these two regions
                    shadeWidthj = 0.03;
                    startFlagj = 0;  %flag=1 if the injury occurred before the recording start, so ensure it is 0 otherwise
                else %otherwise injury or treatment before all regions
                    injuryIndj = locStartBetween/2; %the index of the injury/treatment is between these two regions
                    if injuryIndj == 0.5; injuryIndj = 0.3; end %shift further from violins if injury/treatment before recording start
                    shadeWidthj = 0.03;
                    startFlagj = 1; %flag that the injury occurred before the recording start
                end
            elseif isempty(locStartBetween) %if the treatment/injury was  after the region
                injuryIndj = locEndBetween+ 0.5; %the index of the injury/treatment is between these two regions
                shadeWidthj = 0.03;
                startFlagj = 0; %flag=1 if the injury occurred before the recording start, so ensure it is 0 here
            end

        end %end check for simultaneous regions and injuries/treatments

        if length(injuryIndj)==1 %if there is only one region, or one space regions before which injury occurred..
            InjuryInds(1,j) = injuryIndj; %save the position of injury/treatment j in the list
            InjuryInds(2,j) = shadeWidthj; %save the width of the shading box for injury/treatment j in the list
            InjuryInds(3,j) = startFlagj; %save the width of the shading box for injury/treatment j in the list

        else %otherwise, the injury was administered over multiple regions...
            InjuryInds(1,j) = (injuryIndj(end)-injuryIndj(1))/2+injuryIndj(1); %save the position of injury/treatment j in the list
            numRegionsInvolved = length(injuryIndj); %get the number of regions in the treatment (to determine shading width)
            InjuryInds(2,j) = numRegionsInvolved/2; %save the width of the shading box for injury/treatment j in the list
            InjuryInds(3,j) = startFlagj; %save the width of the shading box for injury/treatment j in the list
        end

    end %end cycle through treatments/injuries

end %end check if there was an injury/treatment


%% Get violin arrays for each region
for rr = 1:AnalysisRegions.numRegions %for each region...

    reg_field = strcat('Region',num2str(rr)); %get the region name to access the data structure
    if plotProps.includeAutocorrelogramsInStatistics == 1 %if user wants to include autocorrelograms in the violins...
        leaderProbArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProb'; %store the metric in the data structure
        numPeaksArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks'; %store the metric in the data structure
        uniformityArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTest'; %store the metric in the data structure
        b1 = zeros(size(RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProbMat));
    else
        b1=eye(size(RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProbMat)); %find the matrix diagonal/autocorrelograms
    end

    b2=ones(ProcessedData.CellCount,ProcessedData.CellCount); b2(:,RecordingMetrics.CorrelogramMetrics.refIndicies)=0;
    b = b1 | b2; %get only the correlograms you want (if specific references were selected or if autocorrelograms were ignored);
 
    %get all cross-correlogram properties (not including autocorrelograms)
    leaderProbArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProbMat(~b); %store the metric in the data structure
    numPeaksArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks(~b); %store the metric in the data structure
    uniformityArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTest(~b); %store the metric in the data structure
    uniformitypArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTestPValue(~b); %store the metric in the data structure
    PeakTimeArray = cell2mat(cellfun(@(x) x, RecordingMetrics.CorrelogramMetrics.(reg_field).correlogramPeakLocations(~b),'UniformOutput',false))'; %list peak times
    PeakLocationArray = cell2mat(cellfun(@(x) 1./abs(x), RecordingMetrics.CorrelogramMetrics.(reg_field).correlogramPeakLocations(~b),'UniformOutput',false))'; %list peak frequencies
    PeakAtZero = isinf(PeakLocationArray); PeakLocationArray = PeakLocationArray(~PeakAtZero); %remove peaks that occurred at t=0, which have f = 1/0 = infinity (which is not physiological)
    PkArray(1:length(PeakLocationArray),rr) = PeakLocationArray; %store peak frequency in the data structure
    PkTArray(1:length(PeakTimeArray),rr) = PeakTimeArray; %store peak times in the data structure
    %double check that there are no entries with all NaNs (for example, if automated region selection was used, and one or more regions had no firing)
    if sum(isnan(leaderProbArray(:,rr))) == ProcessedData.CellCount^2-ProcessedData.CellCount %if the entire matrix/region is empty of events

        %violinplot() is the plotting function used, and can't accept all NaNs.  Therefore, set empty entries to junk/unused values so the violin plots still work, but empty regions are not visible
        leaderProbArray(:,rr) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;
        numPeaksArray(:,rr) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;
        uniformityArray(:,rr) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;
        PkArray(:,rr) = -1000000;
        PkTArray(:,rr) = -1000000; 

    end %end check for empty regions

end %end cycle through regions

PkArray(PkArray==0)=NaN;

%% Plot Violins for each region
figure()
hold on
plot([-3 AnalysisRegions.numRegions+4],[0.5 0.5],'--k','LineWidth',plotProps.lineWidth/2)
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[0 max(leaderProbArray(:))+1],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[0 max(leaderProbArray(:))+1],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[0 max(leaderProbArray(:))+1 max(leaderProbArray(:))+1 0],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,max(leaderProbArray(:))+.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,max(leaderProbArray(:))+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
violinplot(leaderProbArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
hold off
set(gca,'XTick',1:AnalysisRegions.numRegions,'XTickLabel',groupOrder)
set(gca,'YTick',0:.1:1)
ylim([0 1.1])
set(gca,'FontSize',plotProps.FigureFontSize)
box on
grid on
ylabel('Correlogram Area Left of Zero')
xlim([0 AnalysisRegions.numRegions+1])



figure()
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-1 max(numPeaksArray(:))+10],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 max(numPeaksArray(:))+10],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 max(numPeaksArray(:))+10 max(numPeaksArray(:))+10 -1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,max(numPeaksArray(:))+.5,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,max(numPeaksArray(:))+.5,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
violinplot(numPeaksArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
hold off
set(gca,'XTick',1:AnalysisRegions.numRegions,'XTickLabel',groupOrder)
ylim([0 max(numPeaksArray(:))+1])
set(gca,'FontSize',plotProps.FigureFontSize)
if ceil(log10(max(numPeaksArray(:)))) < 2
    set(gca,'YTick',0:max(numPeaksArray(:))+1)
else
    set(gca,'YTick',0:10:max(numPeaksArray(:))+10)
end
box on
grid on
ylabel('Number of Correlogram Peaks')
xlim([0 AnalysisRegions.numRegions+1])


figure()
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-plotProps.correlogramBinMax*1.2 plotProps.correlogramBinMax*1.2],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-plotProps.correlogramBinMax*1.2 plotProps.correlogramBinMax*1.2],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-plotProps.correlogramBinMax*1.2 plotProps.correlogramBinMax*1.2 plotProps.correlogramBinMax*1.2 -plotProps.correlogramBinMax*1.2],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,max(PkTArray(:))*1.1,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,max(PkTArray(:))*1.1,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
violinplot(PkTArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
hold off
set(gca,'XTick',1:AnalysisRegions.numRegions,'XTickLabel',groupOrder)
ylim([-plotProps.correlogramBinMax*1.1 plotProps.correlogramBinMax*1.15])
set(gca,'FontSize',plotProps.FigureFontSize)
set(gca,'YTick',-plotProps.correlogramBinMax:plotProps.correlogramBinMax/10:plotProps.correlogramBinMax)
box on
grid on
ylabel('Correlogram Peak Times (s)')
xlim([0 AnalysisRegions.numRegions+1])


if plotProps.classifyPeakTimesAsFrequencies == 1
    %Get y-axis labels from the min and max possible frequency (based on correlogram bins)
    maxFreq = 1/((Correlograms.CorrelogramBins(end)-Correlograms.CorrelogramBins(1))/length(Correlograms.CorrelogramBins));
    minFreq = 2/(Correlograms.CorrelogramBins(end)-Correlograms.CorrelogramBins(1));
    yTickPts = log10(minFreq/10):1:ceil(1+log10(maxFreq)); yTickPts = 10.^yTickPts;

    figure()
    hold on
    if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
        for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
            plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[minFreq*0.5 10^(1+floor(log10(maxFreq)))],'--k','LineWidth',plotProps.lineWidth/2)
            plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[minFreq*0.5 10^(1+floor(log10(maxFreq)))],'--k','LineWidth',plotProps.lineWidth/2)
            fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[minFreq*0.5 10^(1+floor(log10(maxFreq))) 10^(1+floor(log10(maxFreq))) minFreq*0.5],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
            if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
                text(InjuryInds(1,j) ,maxFreq*1.3,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            else %otherwise injury/treatment is before recording start
                text(InjuryInds(1,j) ,maxFreq*1.3,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            end
        end
    end
    violinplot(PkArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
    hold off
    set(gca,'XTick',1:AnalysisRegions.numRegions,'XTickLabel',groupOrder)
    set(gca,'FontSize',plotProps.FigureFontSize)
    set(gca,'YTick',yTickPts)
    ylim([minFreq*0.9 maxFreq*1.5])
    box on
    grid on
    ylabel('Correlogram Oscillation Frequencies (Hz)')
    xlim([0 AnalysisRegions.numRegions+1])
    set(gca,'YScale','log')
end

figure()
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-1 max(uniformityArray(:))+3],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1  max(uniformityArray(:))+3],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1  max(uniformityArray(:))+3 max(uniformityArray(:))+3 -1 ],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,1.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,1+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
violinplot(uniformityArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
hold off
set(gca,'XTick',1:AnalysisRegions.numRegions,'XTickLabel',groupOrder)
set(gca,'YTick',[0 1],'YTickLabel',{'Nonuniform','Uniform'})
ylim([-0.1 1.1])
set(gca,'FontSize',plotProps.FigureFontSize)
box on
grid on
ylabel('Correlogram Uniformity Test Outcome')
xlim([0 AnalysisRegions.numRegions+1])


figure()
hold on
%plot([-3 AnalysisRegions.numRegions+4],[log10(plotProps.unifPThresh) log10(plotProps.unifPThresh)],'--k','LineWidth',plotProps.lineWidth/2)
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[floor(log10(eps))-1 10],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[floor(log10(eps))-1 10],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[floor(log10(eps))-1 10 10 floor(log10(eps))-1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,1.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,1.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
violinplot(log10(max(eps,uniformitypArray)),1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
plot([-3 AnalysisRegions.numRegions+4],[log10(plotProps.unifPThresh) log10(plotProps.unifPThresh)],'--k','LineWidth',plotProps.lineWidth/2)
hold off
set(gca,'XTick',1:AnalysisRegions.numRegions,'XTickLabel',groupOrder)
set(gca,'YTick',-1000:1:10)
ylim([log10(max(eps,min(uniformitypArray(:))))-1 1.4])
set(gca,'FontSize',plotProps.FigureFontSize)
box on
grid on
ylabel('log_{10}(Probability a Correlogram is Uniform)')
xlim([0 AnalysisRegions.numRegions+1])


%% Get arrays to plot changes between regions
if sw.SelectReferenceCells == 0 && sw.plotDifferencesInCorrelogramMetricsBetweenRegions == 1

    leaderProbDiffArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
    numPeaksDiffArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
    uniformityDiffArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
    for rr = 2:AnalysisRegions.numRegions %for all but the first analysis region...
        for rr2 = 1:rr-1 %cycle through regions before the current region (after will be compared to the current rr the next time rr advances)...
            reg_field = strcat('Region',num2str(rr)); %convert the name of the region to a variable name to call the correct correlogram substructure
            reg_field2 = strcat('Region',num2str(rr2)); %convert the name of the region to a variable name to call the correct correlogram substructure

            diffname = strcat(reg_field,'_Minus_'); diffname = strcat(diffname,reg_field2); %get a variable name for the difference between the two regions

            RecordingMetrics.CorrelogramMetrics.RegionDifferences.(diffname); %save the difference name

            b1=eye(size(RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProbMat)); %find the matrix diagonal/autocorrelograms
            b2=ones(ProcessedData.CellCount,ProcessedData.CellCount);b2(:,RecordingMetrics.CorrelogramMetrics.refIndicies)=0;
            b = b1 | b2; %get only the reference signals you want (if specific references were selected);

            %get all cross-correlogram properties (not including autocorrelograms or specific signals as references if the user has specified so)
            leaderProbDiffArray(:,rr) = RecordingMetrics.CorrelogramMetrics.RegionDifferences.(diffname).leaderProbDifference(~b); %store the metric in the data structure
            numPeaksDiffArray(:,rr) = RecordingMetrics.CorrelogramMetrics.RegionDifferences.(diffname).peakCountDifference(~b); %store the metric in the data structure
            uniformityDiffArray(:,rr) = RecordingMetrics.CorrelogramMetrics.RegionDifferences.(diffname).uniformityDifference(~b); %store the metric in the data structure

            %double check that there are no entries with all NaNs (for example, if automated region selection was used, and one or more regions had no firing)
            if sum(isnan(leaderProbDiffArray(:,rr))) == ProcessedData.CellCount^2-ProcessedData.CellCount %if the entire matrix/region is empty of events

                %violinplot() is the plotting function used, and can't accept all NaNs.  Therefore, set empty entries to junk/unused values so the violin plots still work, but empty regions are not visible
                leaderProbDiffArray(:,rr) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;
                numPeaksDiffArray(:,rr) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;
                uniformityDiffArray(:,rr) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;

            end %end check for empty regions

        end %end cycle through previous regions
    end %end cycle through later regions

    %Pad region 1 with blank placeholder so the x coordinate is correct in the violin plots
    leaderProbDiffArray(:,1) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;
    numPeaksDiffArray(:,1) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;
    uniformityDiffArray(:,1) = zeros(1,ProcessedData.CellCount.^2-ProcessedData.CellCount)-1000000;


    figure()
    hold on
    plot([-3 AnalysisRegions.numRegions+4],[0 0],'--k','LineWidth',plotProps.lineWidth/2)
    if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
        for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
            plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[min(leaderProbDiffArray(:))-1 max(leaderProbDiffArray(:))+1],'--k','LineWidth',plotProps.lineWidth/2)
            plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[min(leaderProbDiffArray(:))-1 max(leaderProbDiffArray(:))+1],'--k','LineWidth',plotProps.lineWidth/2)
            fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[min(leaderProbDiffArray(:))-1 max(leaderProbArray(:))+1 max(leaderProbArray(:))+1 min(leaderProbDiffArray(:))-1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
            if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
                text(InjuryInds(1,j) ,max(leaderProbArray(:))+.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            else %otherwise injury/treatment is before recording start
                text(InjuryInds(1,j) ,max(leaderProbArray(:))+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            end
        end
    end
    violinplot(leaderProbDiffArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
    hold off
    set(gca,'XTick',2:AnalysisRegions.numRegions,'XTickLabel',{groupOrderDelta{2:end}})
    set(gca,'YTick',-1:.1:1)
    ylim([-1.1 1.1])
    set(gca,'FontSize',plotProps.FigureFontSize)
    box on
    grid on
    ylabel('Change in Correlogram Area Left of Zero Between Regions')
    xlim([0 AnalysisRegions.numRegions+1])


    figure()
    hold on
    plot([-3 AnalysisRegions.numRegions+4],[0 0],'--k','LineWidth',plotProps.lineWidth/2)
    if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
        for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
            plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[min(numPeaksDiffArray(:))-1 max(numPeaksDiffArray(:))+3],'--k','LineWidth',plotProps.lineWidth/2)
            plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[min(numPeaksDiffArray(:))-1 max(numPeaksDiffArray(:))+3],'--k','LineWidth',plotProps.lineWidth/2)
            fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[min(numPeaksDiffArray(:))-1 max(numPeaksDiffArray(:))+1 max(numPeaksDiffArray(:))+1 min(numPeaksDiffArray(:))-1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
            if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
                text(InjuryInds(1,j) ,max(numPeaksDiffArray(:))+.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            else %otherwise injury/treatment is before recording start
                text(InjuryInds(1,j) ,max(numPeaksDiffArray(:))+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            end
        end
    end
    violinplot(numPeaksDiffArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
    hold off
    set(gca,'XTick',2:AnalysisRegions.numRegions,'XTickLabel',{groupOrderDelta{2:end}})
    set(gca,'YTick',min(numPeaksDiffArray(abs(numPeaksDiffArray)<1000))-1 : 1: max(numPeaksDiffArray(:))+1)
    ylim([min(numPeaksDiffArray(abs(numPeaksDiffArray)<1000))-1 max(numPeaksDiffArray(:))+1])
    set(gca,'FontSize',plotProps.FigureFontSize)
    box on
    grid on
    ylabel('Change in Correlogram Peak Number Between Regions')
    xlim([0 AnalysisRegions.numRegions+1])


    figure()
    hold on
    plot([-3 AnalysisRegions.numRegions+4],[0 0],'--k','LineWidth',plotProps.lineWidth/2)
    if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
        for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
            plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-2 2],'--k','LineWidth',plotProps.lineWidth/2)
            plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-2 2],'--k','LineWidth',plotProps.lineWidth/2)
            fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-2 2 2 -2],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
            if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
                text(InjuryInds(1,j) ,1.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            else %otherwise injury/treatment is before recording start
                text(InjuryInds(1,j) ,1.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
            end
        end
    end
    violinplot(uniformityDiffArray,1:AnalysisRegions.numRegions,'ViolinColor',{plotProps.regionColors(1:rr,:)});
    hold off
    set(gca,'XTick',2:AnalysisRegions.numRegions,'XTickLabel',{groupOrderDelta{2:end}})
    set(gca,'YTick',-1:1:1,'YTickLabel',{'Lost Uniformity','No Change','Gained Uniformity'})
    ylim([-1.1 1.1])
    set(gca,'FontSize',plotProps.FigureFontSize)
    box on
    grid on
    ylabel('Change in Correlogram Uniformity Between Regions')
    xlim([0 AnalysisRegions.numRegions+1])

end

end %end plotting function