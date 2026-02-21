function [ProcessedData, RecordingMetrics] = calculateRasterMetrics(ProcessedData,plotProps)
%This function calculates and plots recording statistics (raster, firing 
%count histogram, firing interval histogram, and number of raster 
%events/min) independent of correlograms

%Display general information about the recording
disp(' ')
disp(ProcessedData.Comment) %display any comments in the file
disp(['In this recording, ' num2str(ProcessedData.NumberOfActiveElectrodes) ' of ' num2str(ProcessedData.TotalNumberOfElectrodes) ' electrodes had activity and ' num2str(ProcessedData.CellCount) ' unique signals were recorded.'])
disp('Recording Parameters:')
disp(['     Recording Sampling Frequency = ' num2str(ProcessedData.ExperimentSamplFreq) ' Hz'])
if length(unique(ProcessedData.ExperimentTotalGain)) == 1 %if all electrodes had the same gain
    disp(['     Total Gain = ' num2str(ProcessedData.ExperimentTotalGain(1)) ' for all electrodes.'])
else %otherwise, different electrodes had different gain, so display all
    allElectrodez = 1:ProcessedData.TotalNumberOfElectrodes; %list all electrodes (active or otherwise)
    electrodez = unique(cell2mat(cellfun(@(x) str2double(x(2:end-1)), ProcessedData.CellIDs,'UniformOutput',false))); %get each active electrode
    diffGains = unique(ProcessedData.ExperimentTotalGain); %find the different gains on each electrode
    for gz = 1:length(diffGains) %for each gain...
        b = (ProcessedData.ExperimentTotalGain==diffGains(gz)) & ismember(allElectrodez,electrodez); %find electrodes with the given gain
        withGivenGain = allElectrodez(b);
        disp(['     Total Gain = ' num2str(diffGains(gz)) ' for electrode(s): ' num2str(withGivenGain)])
    end
end

%Now obtain overall stats for plotting
FiringCountOfEachCell = cellfun('length',ProcessedData.Raster); %get the total number of times each signal fired
SpikeIntervals = cellfun(@(x) x(2:end)-x(1:end-1),ProcessedData.Raster,'UniformOutput',false); %get the time between spikes

%Save both the overall firing count and the event intervals in the ProcessedData structure
ProcessedData.FiringCountOfEachCell = FiringCountOfEachCell;
ProcessedData.SpikeIntervals = SpikeIntervals;

%% Plot the Rasters
X = cellfun(@(x) x./60,ProcessedData.Raster,'UniformOutput',false); %x axis will be event times
Sizez = cellfun(@size,ProcessedData.Raster,'UniformOutput',false); %get the length of the raster for each signal
for i = 1:ProcessedData.CellCount %cycle through signals...
    IDArr{i} = i+zeros(Sizez{i}); %make y coordinate the signal ID and make as long as the signal's raster
end %end cycle through signals

figure()
hold on
p = cellfun(@plot,X,IDArr); %plot rasters (y axis is signal, x axis is time)
for i = 1:ProcessedData.CellCount %cycle through signals...
    p(i).Marker = '|';
    p(i).LineStyle = 'none';
    p(i).Color = [0 0 0];
end %end set raster mark appearance
hold off
set(gca,'YTick',1:ProcessedData.CellCount,'YTickLabel',ProcessedData.CellIDs)
xlim([0 ProcessedData.ExperimentEndTime./60]) %convert time (s) to min
ylim([0 ProcessedData.CellCount+1])
ylabel('Signal')
xlabel('Time (min)')
title('Signal Rasters')
box on
set(gca,'FontSize',plotProps.FigureFontSize)


%% Plot Raster-Related Histograms
bins = linspace(0,1+log10(max(FiringCountOfEachCell)),round(ProcessedData.CellCount/1.5)); bins = 10.^bins;
figure()
histogram(FiringCountOfEachCell,bins,'FaceAlpha',0.3)
xlabel('Number of Times a Signal Fired During the Recording')
ylabel('Number of Signals')
title('Firing Count Histogram')
box on
grid on
set(gca,'XScale','log')
set(gca,'FontSize',plotProps.FigureFontSize)
set(gca,'XTick',10.^[0:1:round(log10(max(bins)))])
xlim(10.^[0 round(log10(max(bins)))])

bins = linspace(-5,4,100); bins = 10.^bins;
figure()
histogram(vertcat(SpikeIntervals{:}),bins,'FaceAlpha',0.3)
xlabel('Interval Between Events (s)')
ylabel('Count')
title('Event Interval Histogram')
set(gca,'XScale','log')
box on
grid on
set(gca,'FontSize',plotProps.FigureFontSize)
set(gca,'XTick',10.^[round(log10(min(bins))-1):1:round(log10(max(bins)))])
xlim(10.^[round(log10(min(bins))-1) round(log10(max(bins)))])

%% Save spikes/min vs. time

%Make bins of time for each minute
numMinutes = ProcessedData.ExperimentEndTime/60; %convert end time from s to min
minBins = 0:numMinutes;

spikesPerMin = cellfun(@(x) histcounts( x./60, minBins),ProcessedData.Raster,'UniformOutput',false); %get the spikes/minute for each signal

meanSPM = mean(vertcat(spikesPerMin{:}),1); %calculate mean spikes/minute for the entire network
stdSPM = std(vertcat(spikesPerMin{:}),1); %calculate standard deviation in spikes/min for the entire network


%save the spikes per minute metrics
RecordingMetrics.meanSPM = meanSPM;
RecordingMetrics.stdSPM = stdSPM;
RecordingMetrics.minBinsSPM = minBins(2:end);
RecordingMetrics.individualSPM = spikesPerMin;
RecordingMetrics.minBins = minBins;

end %end function
