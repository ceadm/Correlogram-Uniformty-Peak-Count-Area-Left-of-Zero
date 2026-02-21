function summarizeCorrelograms(sw, plotProps, AnalysisRegions, InjuryIndicies, ProcessedData, RecordingMetrics, Correlograms)
%This function classifies correlograms based on metrics (uniformity, peak 
%count, area left of zero) and plots classification-related outputs

%% Find where recording regions occurr relative to injuries or treatments
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury/treatment in the history of the recording

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

        if any(simultaneousRegAndInj) %if any row had both positive and negative times/occurred during an injury
            injuryIndj = find(simultaneousRegAndInj);
            shadeWidthj = 0.4;
            startFlagj = 0; %flag=1 if the injury occurred before the recording start

        else %otherwise the regions are either completely after or completely after the given injury/treatment
            locStartBetween = find(endNegSum==0,1,'first'); %find the first region after the injury/treatment
            locEndBetween = find(endPosSum==0,1,'last'); %find the last region before the injury/treatment

            if ~isempty(locEndBetween) && ~isempty(locStartBetween) %if the injury/treatment was during the recording
                injuryIndj = locStartBetween-(locStartBetween - locEndBetween)/2; %the index of the injury/treatment is between these two regions
                shadeWidthj = 0.03;
                startFlagj = 0;  %flag=1 if the injury occurred before the recording start

            elseif isempty(locEndBetween) %if the treatment/injury was before the recording started or after all selected regions
                if InjuryIndicies.InjuryStartMin(j) > 0 %injury was after all analyzed regions
                    injuryIndj = locStartBetween-.5; %the index of the injury/treatment is between these two regions
                    shadeWidthj = 0.03;
                    startFlagj = 0;  %flag=1 if the injury occurred before the recording startstartFlagj = 0; %flag=1 if the injury occurred before the recording start
                else %otherwise injury or treatment before all regions
                    injuryIndj = locStartBetween/2; %the index of the injury/treatment is between these two regions
                    if injuryIndj == 0.5; injuryIndj = 0.3; end %shift further from violins if injury/treatment before recording start
                    shadeWidthj = 0.03;
                    startFlagj = 1; %flag that the injury occurred before the recording start
                end
            elseif isempty(locStartBetween) %if the treatment/injury was  after the region
                injuryIndj = locEndBetween+ 0.5; %the index of the injury/treatment is between these two regions
                shadeWidthj = 0.03;
                startFlagj = 0; %flag=1 if the injury occurred before the recording start
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

%% Get recording zones
%Zones = times relative to injury or treatment (i.e. before treatment, between treatments, after treatment)
groupOrder = ProcessedData.RegionStats.groupOrder;
timeRange = ProcessedData.RegionStats.timeRange;

% Sort the analysis regions into each recording zone
Zonez = zeros(1,AnalysisRegions.numRegions); %row = zone, column = region
for zz = 1:size(timeRange,2) %cycle through zones...

    startTimeLimit = timeRange(1,zz); %the zone starts mid recording, so just look for anything from zone start to zone end
    endTimeLimit = timeRange(2,zz); %the zone ends mid recording as well

    %Find analysis regions within the zone
    %Note that any regions that are part of multiple zones are not marked as being part of either zone
    %however, such regions will be used as comparisons against those in the zone later

    afterZoneStart = AnalysisRegions.Seconds>=startTimeLimit; %find analysis region times after zone start
    afterZoneStart = sum(afterZoneStart,2)==2; %get analysis regions where the region start and end are both after the start of the zone

    beforeZoneEnd = AnalysisRegions.Seconds<=endTimeLimit; %find analysis region times before zone end
    beforeZoneEnd = sum(beforeZoneEnd,2)==2; %get analysis regions where the region start and end are both before the end of the zone

    inZoneInds = afterZoneStart&beforeZoneEnd; %find regions contained within the zone

    if ~any(inZoneInds==1) %if the zone is smaller than a single analysis region, ensure that those regions (where the treatmment is administered) are added to the zone
        disp(['     Zone ' num2str(zz) ' (' char(groupOrder{zz}) ') began and ended in a single region or in the middle of two adjacent regions.'])

        afterZoneStart = AnalysisRegions.Seconds>=startTimeLimit; %find analysis region times after zone start
        afterZoneStart = sum(afterZoneStart,2)==1; %get analysis regions where the region start and end are both after the start of the zone

        beforeZoneEnd = AnalysisRegions.Seconds<=endTimeLimit; %find analysis region times before zone end
        beforeZoneEnd = sum(beforeZoneEnd,2)==1; %get analysis regions where the region start and end are both before the end of the zone

        zoneLimitz = find(afterZoneStart):find(beforeZoneEnd);
        inZoneInds = false(size(afterZoneStart)); inZoneInds(zoneLimitz) = true; %find regions contained within the zone
    end
    Zonez(inZoneInds) = zz; %mark which zone the region belongs to

end %end cycle through zones

numZones = length(unique(Zonez(Zonez>0))); %count the number of recording zones



%% Get correlogram data arrays

%Determine whether to exclude autocorrelograms from the analysis based on user preferences
if plotProps.includeAutocorrelogramsInStatistics == 1 %if user wants to include autocorrelograms in the violins...
    b1 = zeros(size(RecordingMetrics.CorrelogramMetrics.Region1.leaderProbMat)); %include all indicies
else %otherwise, only include crosscorrelograms in violins
    %get the autocorrelograms (same number as the signal count, unless only select references are specified)
    b1=eye(size(RecordingMetrics.CorrelogramMetrics.Region1.leaderProbMat)); %find the matrix diagonal/autocorrelograms
end

%If only selecting a subset of signals as the reference signal, take only the relevent indicies (of those comparisons only)
b2=ones(ProcessedData.CellCount,ProcessedData.CellCount); b2(:,RecordingMetrics.CorrelogramMetrics.refIndicies)=0;

b=(b1 & ~b2); autoCorNum = sum(b(:)); %count how many correlograms (if any) are being excluded

%Initialize arrays
leaderProbArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
numPeaksArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
uniformityArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);
uniformitypArray = zeros(ProcessedData.CellCount.*length(RecordingMetrics.CorrelogramMetrics.refIndicies)-autoCorNum,AnalysisRegions.numRegions);

%For peak positions (which provide frequencies for any rhythmic activity), each region can have a different number of peaks.  Therefore, regions with less than the maximum number of peaks need to be padded with dummy numbers so that they are all the same size for violin plot
MaxNumPk = 0;
for rr = 1:AnalysisRegions.numRegions %find the maximum possible number of peaks for all regions
    reg_field = strcat('Region',num2str(rr)); %convert the name of the region to a variable name to call the correct correlogram substructure
    MaxNumPk = max(MaxNumPk,max(RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks));
end
PkArray = NaN(MaxNumPk,AnalysisRegions.numRegions);  %initialize peak frequency array

SourceCorrelogram = zeros(size(PkArray)); %initialize array to store an ID for the correlogram in which each peak is located

%% Populate the correlogram arrays for each analysis region
for rr = 1:AnalysisRegions.numRegions %for each region...

    reg_field = strcat('Region',num2str(rr)); %get the region name to access the data structure
    if plotProps.includeAutocorrelogramsInStatistics == 1 %if user wants to include autocorrelograms in the violins
        leaderProbArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProb'; %store the metric in the data structure
        numPeaksArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks'; %store the metric in the data structure
        uniformityArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTest'; %store the metric in the data structure
        b1 = zeros(size(RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProbMat));
    else
        b1=eye(size(RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProbMat)); %find the matrix diagonal/autocorrelograms
    end

    b2=ones(ProcessedData.CellCount,ProcessedData.CellCount); b2(:,RecordingMetrics.CorrelogramMetrics.refIndicies)=0;
    b = b1 | b2; %get only the reference signals you want (if specific references were selected);

    %get all cross-correlogram properties (not including autocorrelograms or excluded references if specified)
    leaderProbArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProbMat(~b); %store the metric in the data structure
    numPeaksArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks(~b); %store the metric in the data structure
    uniformityArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTest(~b); %store the metric in the data structure
    uniformitypArray(:,rr) = RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTestPValue(~b); %store the metric in the data structure
    
    if plotProps.classifyPeakTimesAsFrequencies == 1 %if classifying peak times according to frequency (Hz, f = 1/t)...
        PeakLocationArray = cell2mat(cellfun(@(x) 1./abs(x), RecordingMetrics.CorrelogramMetrics.(reg_field).correlogramPeakLocations(~b),'UniformOutput',false))'; %list peak frequencies
        PeakAtZero = isinf(PeakLocationArray); PeakLocationArray = PeakLocationArray(~PeakAtZero); %remove peaks that occurred at t=0, which have f = 1/0 = infinity (which is not physiological)
        PkArray(1:length(PeakLocationArray),rr) = PeakLocationArray; %store the metric in the data structure
    else %otherwise classify peak times accorsing to time (s)
        PeakLocationArray = cell2mat(cellfun(@(x) abs(x), RecordingMetrics.CorrelogramMetrics.(reg_field).correlogramPeakLocations(~b),'UniformOutput',false))'; %list peak frequencies
        PkArray(1:length(PeakLocationArray),rr) = PeakLocationArray; %store the metric in the data structure
    end

    %Store the ID of the correlogram from which each peak originated
    startInd = 1; %start index for source correlogram IDs (since the number of peaks can be > 1 for any correlogram, the peak frequency array is larger than the correlogram arrays.  Therefore, it is necessary to track which correlogram contained each peak)
    for i = 1:ProcessedData.CellCount^2 %cycle through correlograms...
        N = RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks(i); %Get the number of correlogram peaks
        N2 = sum(isinf(1./(abs(RecordingMetrics.CorrelogramMetrics.(reg_field).correlogramPeakLocations{i})))); %determine if f = Inf (unphysiological) for any of the peaks

        if isnan(N) || N > N2 %if there is at least one peak not equal to infinity...
            if isnan(N) %set the value of N if the correlogram is empty
                N = 1; N2 = 0;
            end

            SourceCorrelogram(startInd:startInd+(N-N2)-1,rr) = i; %save correlogram source ID in the matrix
            startInd = startInd+(N-N2); %advance to the start for the next correlogram ID

        end %end check that there is at least one peak not equal to infinity

    end %end cycle through correlograms

end %end cycle through regions

PkArray(PkArray==0)=NaN;


%% Classify leader/follower nature of the correlograms
leaderProbGroupBoundaries = 0:.2:1;
leaderClassificationNames = {'Weak', 'Fairly Weak', 'Intermediate', 'Fairly Strong', 'Strong'};
NumberOfleaderClassificationGroups = length(leaderClassificationNames); %get all possible groups for nonempty correlograms

fracInProbGroup = NaN(size(leaderProbGroupBoundaries,2)-1,AnalysisRegions.numRegions); %Row = grouping, Col = Recording Region
LFGroupAssignments = NaN(size(leaderProbArray)); %initialize array to store leader/follower classifications

StrongWeakClass = abs(leaderProbArray - 0.5)./0.5; %Bin into five categories of weak->strong leader/follower characteristic.  A value of 1 is strong, 0 is weak

for g = 1:length(leaderProbGroupBoundaries)-1
    %Get bounds on the probabilities in range
    a = leaderProbGroupBoundaries(g); b = leaderProbGroupBoundaries(g+1);

    %Find the values of how strong/week a given correlogram leads/follows within in the given range (a - b)
    if g == 1 %at the start of category listing, so pad the zero (lower group bound) to avoid machine/rounding error
        ProbInBin = (StrongWeakClass>=a-1 & StrongWeakClass<b); %find leader/follower strength classifications within the given range
    elseif b~=1 %as long as b is not the end of the list
        ProbInBin = (StrongWeakClass>=a & StrongWeakClass<b); %find leader/follower strength classifications within the given range
    else %otherwise, end of the category listing, so be sure to include the highest limit (don't reserve for the next group as with non-end elements of the categories)
        %at the start of category listing, so pad the 1 (upper group bound) to avoid machine/rounding error
        ProbInBin = (StrongWeakClass>=a & StrongWeakClass<=b+1); %find leader/follower strength classifications within the given range
    end

    LFGroupAssignments(ProbInBin) = g; %save the logical array marking which correlograms fall into which bin

    fracInProbGroup(g,1:end) = sum(ProbInBin,1); %sum up how many correlograms fit in the given bin

end

%Do not count correlograms with 0 events in the total
EmptyCorrelograms = isnan(leaderProbArray); %Find correlograms with 0 events
totalNumberOfCorrelograms = sum(~EmptyCorrelograms,1); %Get the actual/non-empty count

fracInProbGroup = fracInProbGroup./totalNumberOfCorrelograms; %Get the fraction of correlograms in each bin


%% Classify peak counts of the correlograms
pkCountGroupBoundaries = min(numPeaksArray(:)):1:min(plotProps.maxPeakCountBeforeNoise,max(numPeaksArray(:))); %bin peak counts from the minimum to a user-specified maximum or the true maximum (whichever max value is smaller)
pkCountNames = string(pkCountGroupBoundaries); %group names are just the peak count

if plotProps.maxPeakCountBeforeNoise < max(numPeaksArray(:)) %the last group/peak count is a user-specified maximum, which may be smaller than the true maximum, so add "or more" to the name if this is the case
    pkCountNames{end} = strcat(pkCountNames{end},' or More');
end
NumberOfpkCountClassificationGroups = length(pkCountNames); %get all possible groups for nonempty correlograms

fracInPeakGroup = NaN(size(pkCountGroupBoundaries,2)-1,AnalysisRegions.numRegions); %Row = grouping, Col = Recording Region
PCGroupAssignments = NaN(size(numPeaksArray)); %initialize array to store peak count classifications

for g = 1:length(pkCountGroupBoundaries)
    %Get bounds on the probabilities in range
    if g == length(pkCountGroupBoundaries)
        ProbInBin = numPeaksArray >= pkCountGroupBoundaries(g); %find correlograms with the given number of peaks or more if at the end of the group list
    else
        ProbInBin = numPeaksArray == pkCountGroupBoundaries(g); %find correlograms with the given number of peaks
    end

    PCGroupAssignments(ProbInBin) = g; %save the logical array marking which correlograms fall into which bin

    fracInPeakGroup(g,1:end) = sum(ProbInBin,1); %sum up how many correlograms fit in the given bin

end

%Do not count correlograms with 0 events in the total
EmptyCorrelograms = isnan(numPeaksArray); %Find correlograms with 0 events
totalNumberOfCorrelograms = sum(~EmptyCorrelograms,1); %Get the actual/non-empty count

fracInPeakGroup = fracInPeakGroup./totalNumberOfCorrelograms; %Get the fraction of correlograms in each bin

%% Uniformity is already classified as uniform and nonuniform, so just get the fraction of uniform and nonuniform in each data analysis region
unifGroupBoundaries = [0 1];
unifNames = {'Nonuniform','Uniform'};
NumberOfunifClassificationGroups = length(unifGroupBoundaries); %get all possible groups for nonempty correlograms

fracInUnifGroup = NaN(size(unifGroupBoundaries,2)-1,AnalysisRegions.numRegions); %Row = grouping, Col = Recording Region
unifGroupAssignments = NaN(size(uniformityArray)); %initialize array to store peak count classifications

for g = 1:length(unifGroupBoundaries)
    %Get bounds on the probabilities in range
    ProbInBin = uniformityArray == unifGroupBoundaries(g); %find correlograms with the given uniformity classification

    unifGroupAssignments(ProbInBin) = g; %save the logical array marking which correlograms fall into which bin

    fracInUnifGroup(g,1:end) = sum(ProbInBin,1); %sum up how many correlograms fit in the given bin

end

%Do not count correlograms with 0 events in the total
EmptyCorrelograms = isnan(uniformityArray); %Find correlograms with 0 events
totalNumberOfCorrelograms = sum(~EmptyCorrelograms,1); %Get the actual/non-empty count

fracInUnifGroup = fracInUnifGroup./totalNumberOfCorrelograms; %Get the fraction of correlograms in each bin



%% Classify correlogram peak times or peak frequencies
% EEG Categories (if classifying by frequency)
%     Delta: f < 4
%     Theta: 4 <= f <= 7
%     Alpha: 8 <= f <= 12
%     Beta: 13 <= f <= 30
%     Gamma: f > 32

if plotProps.classifyPeakTimesAsFrequencies == 1 %if classifying peak times according to frequency (Hz, f = 1/t)...
    startF = [0 4 7 8 12 13 30 32]; %list the minimum frequency of each category
    endF = [4 7 8 12 13 30 32 max(PkArray(:))+10]; %list the maximum frequency of each category
    freqClassificationNames = {'Delta Band (<4 Hz)', 'Theta Band (4-7 Hz)', 'Between Theta and Alpha (7-8 Hz)', 'Alpha Band (8-12 Hz)', 'Between Alpha and Beta (12-13 Hz)', 'Beta Band (13-30 Hz)', 'Between Beta and Gamma (30-32 Hz)', 'Gamma Band (>32 Hz)'};
    NumberOfFreqClassificationGroups = length(freqClassificationNames); %count the number of frequency categories
else
    stpowr = [-Inf floor(log10(Correlograms.CorrelogramBins(end)-Correlograms.CorrelogramBins(end-1))):1:ceil(log10(max(PkArray(:))))-1]; %get the log10 power of time bin starts
    startF = 10.^stpowr; %list the minimum time of each category
    endpowr = [floor(log10(Correlograms.CorrelogramBins(end)-Correlograms.CorrelogramBins(end-1))):1:ceil(log10(max(PkArray(:))))]; %get the log10 power of time bin ends
    endF = 10.^endpowr; %list the maximum time of each category

    freqClassificationNames = repmat({' '},1,length(stpowr)); %initialize array to save time grouping names
    for i = 1:length(stpowr)
        if i == 1; ind = i+1; else; ind = i; end %get index to access the start time name

        %get the start and end time strings
        if stpowr(ind) == -12; sttstr = '1 picosecond'; elseif stpowr(ind) == -11; sttstr = '10 picoseconds'; elseif stpowr(ind) == -10; sttstr = '100 picoseconds'; elseif stpowr(ind) == -9; sttstr = '1 nanosecond'; elseif stpowr(ind) == -8; sttstr = '10 nanoseconds'; elseif stpowr(ind) == -7; sttstr = '100 nanoseconds'; elseif stpowr(ind) == -6; sttstr = '1 microsecond'; elseif stpowr(ind) == -5; sttstr = '10 microseconds'; elseif stpowr(ind) == -4; sttstr = '100 microseconds'; elseif stpowr(ind) == -3; sttstr = '1 millisecond'; elseif stpowr(ind) == -2; sttstr = '10 milliseconds'; elseif stpowr(ind) == -1; sttstr = '100 milliseconds'; elseif stpowr(ind) == 0; sttstr = '1 second'; elseif stpowr(ind) == 1; sttstr = '10 seconds'; elseif stpowr(ind) == 2; sttstr = '100 seconds'; elseif stpowr(ind) == 3; sttstr = '1,000 seconds';  elseif stpowr(ind) == 4; sttstr = '10,000 seconds'; elseif stpowr(ind) == 5; sttstr = '100,000 seconds'; elseif stpowr(ind) == 6; sttstr = '1,000,000 seconds'; end
        if endpowr(i) == -12;  endtstr = '1 picosecond'; elseif endpowr(i) == -11; endtstr = '10 picoseconds'; elseif endpowr(i) == -10; endtstr = '100 picoseconds'; elseif endpowr(i) == -9; endtstr = '1 nanosecond'; elseif endpowr(i) == -8;  endtstr = '10 nanoseconds'; elseif endpowr(i) == -7; endtstr = '100 nanoseconds'; elseif endpowr(i) == -6;  endtstr = '1 microsecond'; elseif endpowr(i) == -5;  endtstr = '10 microseconds'; elseif endpowr(i) == -4;  endtstr = '100 microseconds'; elseif endpowr(i) == -3; endtstr = '1 millisecond'; elseif endpowr(i) == -2; endtstr = '10 milliseconds'; elseif endpowr(i) == -1; endtstr = '100 milliseconds'; elseif endpowr(i) == 0; endtstr = '1 second'; elseif endpowr(i) == 1;  endtstr = '10 seconds'; elseif endpowr(i) == 2; endtstr = '100 seconds'; elseif endpowr(i) == 3; endtstr = '1,000 seconds'; elseif endpowr(i) == 4; endtstr = '10,000 seconds'; elseif endpowr(i) == 5; endtstr = '100,000 seconds';  elseif endpowr(i) == 6; endtstr = '1,000,000 seconds';  end
    
        if i == 1 %get the first entry in the categories
            tstr = [sttstr ' or less'];
        else %now get the non-start time ranges
            startchar = strfind(sttstr,' '); startchar = startchar(1); %get the index when the first number stops
            if startchar == 4 %if the first number is 100, use the full start name
                tstr = [sttstr ' - ' endtstr];
            else %otherwise, you only need to write the units once
                tstr = [sttstr(1:startchar-1) ' - ' endtstr];
            end
        end %end get the category range
        freqClassificationNames{i} = tstr; %save the category in the time classification name array
    end

    NumberOfFreqClassificationGroups = length(freqClassificationNames); %count the number of peak time categories
    
end

%Initialize peak frequency (or time) classification arrays
fracInFreqGroup = NaN(NumberOfFreqClassificationGroups,AnalysisRegions.numRegions); %Row = fraction of correlogram peaks in grouping, Col = Recording Region
meanFreqInFreqGroup = NaN(NumberOfFreqClassificationGroups,AnalysisRegions.numRegions); %Row = mean frequency within the given group, Col = Recording Region
sdFreqInFreqGroup = NaN(NumberOfFreqClassificationGroups,AnalysisRegions.numRegions); %Row = standard deviation in the frequency within the given group, Col = Recording Region
seFreqInFreqGroup = NaN(NumberOfFreqClassificationGroups,AnalysisRegions.numRegions); %Row = standard deviation in the frequency within the given group, Col = Recording Region
minFreqInFreqGroup = NaN(NumberOfFreqClassificationGroups,AnalysisRegions.numRegions); %Row = min frequency within the given group, Col = Recording Region
maxFreqInFreqGroup = NaN(NumberOfFreqClassificationGroups,AnalysisRegions.numRegions); %Row = max frequency within the given group, Col = Recording Region
FreqGroupAssignments = NaN(size(PkArray)); %initialize array to store peak frequency classifications

for g = 1:NumberOfFreqClassificationGroups %for each classification grouping...
    
    if g == 1
        b = PkArray<= endF(g);
    else
        b = PkArray>startF(g) & PkArray <= endF(g); %find frequencies that fall within the given classification/grouping
    end
    FreqGroupAssignments(b) = g; %save the logical array marking peaks in the given classification

    fracInFreqGroup(g,:) = sum(b,1); %count how many correlogram peaks have a frequency in the given classification

    %min and max functions will return an empty array if no peaks are in group g for a given region.  This interferes with MATLAB's indexing, so empty regions must be found and ignored for the min and max functions
    regionsWithPeaksInGroup = sum(b,1)>0; %make a logical indicating regions with at least one peak in the given classification group
    regionsWithGroup = 1:AnalysisRegions.numRegions; regionsWithGroup = regionsWithGroup(regionsWithPeaksInGroup); %list regions with at least one peak in the given classification group
    
    if any(regionsWithPeaksInGroup)
        meanFreqInFreqGroup(g,:) = cell2mat(arrayfun(@(col) mean(PkArray(b(:,col),col),'omitnan'),1:AnalysisRegions.numRegions,'UniformOutput',false));  % mean frequency of correlogram peaks in the given classification
        medianFreqInFreqGroup(g,:) = cell2mat(arrayfun(@(col) nanmedian(PkArray(b(:,col),col)),1:AnalysisRegions.numRegions,'UniformOutput',false)); % median frequency of correlogram peaks in the given classification
        sdFreqInFreqGroup(g,:) = cell2mat(arrayfun(@(col) std(PkArray(b(:,col),col),'omitnan'),1:AnalysisRegions.numRegions,'UniformOutput',false));  % standard deviation in the frequencies of correlogram peaks in the given classification
        seFreqInFreqGroup(g,:) = cell2mat(arrayfun(@(col) std(PkArray(b(:,col),col)/sqrt(length(b(:,col))),'omitnan'),1:AnalysisRegions.numRegions,'UniformOutput',false));  % standard error in the frequencies of correlogram peaks in the given classification
        minFreqInFreqGroup(g,regionsWithGroup) = cell2mat(arrayfun(@(col) min(PkArray(b(:,col),col),[],'omitnan'),regionsWithGroup,'UniformOutput',false));  % min frequency of correlogram peaks in the given classification
        maxFreqInFreqGroup(g,regionsWithGroup) = cell2mat(arrayfun(@(col) max(PkArray(b(:,col),col),[],'omitnan'),regionsWithGroup,'UniformOutput',false)); % max frequency of correlogram peaks in the given classification
    end

end

%Do not count correlograms with 0 events in the total for the % of peaks in each classification
EmptyCorrelograms = isnan(PkArray) | PkArray == 0; %Find correlograms with 0 events
totalNumberOfCorrelograms = sum(~EmptyCorrelograms,1); %Get the actual/non-empty count

fracInFreqGroup = fracInFreqGroup./totalNumberOfCorrelograms; %Get the fraction of correlograms in each bin (sum of all rows for each column = 1)


%% Make a structure containing indicies for each classification group
%(for ease of accessing later/to prevent constantly using "find")

for g = 1:NumberOfleaderClassificationGroups
    g_name = strcat('G',num2str(g));
    FindGroups.LeaderFollower.(g_name) =  LFGroupAssignments==g;
end
for g = 1:NumberOfpkCountClassificationGroups
    g_name = strcat('G',num2str(g));
    FindGroups.PeakCount.(g_name) =  PCGroupAssignments==g;
end
for g = 1:NumberOfunifClassificationGroups
    g_name = strcat('G',num2str(g));
    FindGroups.Unif.(g_name) =  unifGroupAssignments==g;
end
for g = 1:NumberOfFreqClassificationGroups
    g_name = strcat('G',num2str(g));
    FindGroups.FreqClassification.(g_name) =  FreqGroupAssignments==g;
end

%% Analyze classification changes between treatments/injuries
%Get the probability that the correlogram now has its current classification, given the classification in the previous recording region
for zz = 1: numZones %cycle through recording zones
    zoneName = strcat('Zone',num2str(zz)); %get the zone name to access the data structure
    startR = find(Zonez==zz,1,'first'); %get the first region in the zone
    endR = find(Zonez==zz,1,'last'); %get the last region in the zone

    LFGroupz = zeros(NumberOfleaderClassificationGroups,NumberOfleaderClassificationGroups); %Initialize array to save probability of being in a group given the past group
    PCGroupz = zeros(NumberOfpkCountClassificationGroups,NumberOfpkCountClassificationGroups); %Initialize array to save probability of being in a group given the past group
    UnifGroupz = zeros(NumberOfunifClassificationGroups,NumberOfunifClassificationGroups); %Initialize array to save probability of being in a group given the past group

    %row = previous group, column = probability of current group
    for rr = max(2,startR):endR %cycle through regions within the zone...

        %Get the probability of a current leader/follower classification given the value of the previous classification
        for g = 1:NumberOfleaderClassificationGroups %cycle through past group...
            g_name = ['G'  num2str(g)]; %get the group name
            indArr = FindGroups.LeaderFollower.(g_name); %find correlograms classified in the given group

            previouslyInGroup = indArr(:,rr-1); %find correlograms belonging to the given group in the previous analysis region
            currentGroup = LFGroupAssignments(previouslyInGroup,rr); %find how those correlograms are classified in the current analysis region

            for g2 = 1:NumberOfleaderClassificationGroups %for each current group...
                LFGroupz(g,g2) = LFGroupz(g,g2)+sum(currentGroup==g2); %add current groups to histogram counts
            end %end cycle through current groups

        end %end cycle through past group


        %Get the probability of a current peak count given the value of the previous peak count
        for g = 1:NumberOfpkCountClassificationGroups %cycle through past group...
            g_name = ['G'  num2str(g)]; %get the group name
            indArr = FindGroups.PeakCount.(g_name); %find correlograms classified in the given group

            previouslyInGroup = indArr(:,rr-1); %find correlograms belonging to the given group in the previous analysis region
            currentGroup = PCGroupAssignments(previouslyInGroup,rr); %find how those correlograms are classified in the current analysis region

            for g2 = 1:NumberOfpkCountClassificationGroups %for each current group...
                PCGroupz(g,g2) = PCGroupz(g,g2)+sum(currentGroup==g2); %add current groups to histogram counts
            end %end cycle through current groups

        end %end cycle through past group


        %Get the probability of a current uniformity classification given the value of the previous uniformity classification
        for g = 1:NumberOfunifClassificationGroups %cycle through past group...
            g_name = ['G'  num2str(g)]; %get the group name
            indArr = FindGroups.Unif.(g_name); %find correlograms classified in the given group

            previouslyInGroup = indArr(:,rr-1); %find correlograms belonging to the given group in the previous analysis region
            currentGroup = unifGroupAssignments(previouslyInGroup,rr); %find how those correlograms are classified in the current analysis region

            for g2 = 1:NumberOfunifClassificationGroups %for each current group...
                UnifGroupz(g,g2) = UnifGroupz(g,g2)+sum(currentGroup==g2); %add current groups to histogram counts
            end %end cycle through current groups

        end %end cycle through past group

    end %end cycle through analysis regions within the zone

    ZoneHistograms.LeaderFollower.(zoneName) = LFGroupz./sum(LFGroupz,2); %normalize histogram to probabilities and save in the structure so you can plot histograms for each zone
    ZoneHistograms.PeakCount.(zoneName) = PCGroupz./sum(PCGroupz,2); %normalize histogram to probabilities and save in the structure so you can plot histograms for each zone
    ZoneHistograms.Unif.(zoneName) = UnifGroupz./sum(UnifGroupz,2); %normalize histogram to probabilities and save in the structure so you can plot histograms for each zone

end %end cycle through recording zones


%% Make Plots of Classification Outputs

%Initialize plot markers for ease of plotting
markerNmArr = {'o','*','x','square','diamond','^','v','>','<','pentagram','hexagram','+','.'}; %Make an array to change the plot marker shapes for each category of classified data
markerSzArr = [plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize*0.4 plotProps.markerSize];

%Make sure there are enough marker indicators for all plots to be made
maxClassificationCount = (1+min(plotProps.maxPeakCountBeforeNoise,max(numPeaksArray(:)))); %The highest number of classifications will be peak counts, so make sure length(markerArray) is larger than the max possible number of peak counts
maxClassificationCount = max(length(ProcessedData.RegionStats.groupOrder),maxClassificationCount); %the highest possible number of lines on that plots will be the number of analysis regions

F = maxClassificationCount/length(markerNmArr);  %The highest number of classifications will be peak counts, so make sure length(markerArray) is larger than the max possible number of peak counts
if F > 1 %if the marker array is too short to plot all peak classifications
    F = ceil(F); %round up so the array will eventually be longer than needed
    for i = 1:F %as many times as needed...
        markerNmArr = [markerNmArr markerNmArr];  %lengthen the marker name array as many times as needed
        markerSzArr = [markerSzArr markerSzArr];
    end
end %end check marker specification is long enough



%Plot the fraction of correlogram peaks in each leader/follower classification
cc = pink(NumberOfleaderClassificationGroups+4); %set plot colors

figure()
%Indicate where injuries/treatments occurred
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2 2 -1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,max(fracInProbGroup(:))+.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,max(fracInProbGroup(:))+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
%Now plot the fraction of correlogram peaks in each classification
h=plot(1:AnalysisRegions.numRegions,fracInProbGroup','.-','LineWidth',plotProps.lineWidth);
hold off
for c = 1:size(h,1)
    h(c).Color = cc(c,:);
    h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
end
set(gca,'XTick',1:AnalysisRegions.numRegions)
if isfield(AnalysisRegions,'RegionLabels')
    set(gca,'XTickLabel',AnalysisRegions.RegionLabels)
end

xlim([0.5 AnalysisRegions.numRegions+0.5])
xlabel('Recording Region')
ylim([0 1.2*(max(fracInProbGroup(:)))])
ylabel('Fraction of Correlograms with a given Classification')
legend(h,leaderClassificationNames)
set(gca,'FontSize',plotProps.FigureFontSize)
title('Classification of Leader/Follower Nature Through Time')
box on
set(gca,'YTick',0:0.05:1)
ax = gca;
ax.YGrid = 'on';



%Plot the fraction of correlogram peaks in each peak count
cc = pink(NumberOfpkCountClassificationGroups+4); %set plot colors

figure()
%Indicate where injuries/treatments occurred
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2 2 -1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,max(fracInPeakGroup(:))+.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,max(fracInPeakGroup(:))+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
%Now plot the fraction of correlogram peaks in each classification
h=plot(1:AnalysisRegions.numRegions,fracInPeakGroup','.-','LineWidth',plotProps.lineWidth);
hold off
for c = 1:size(h,1)
    h(c).Color = cc(c,:);
    h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
end
set(gca,'XTick',1:AnalysisRegions.numRegions)
if isfield(AnalysisRegions,'RegionLabels')
    set(gca,'XTickLabel',AnalysisRegions.RegionLabels)
end

xlim([0.5 AnalysisRegions.numRegions+0.5])
xlabel('Recording Region')
ylim([0 1.2*(max(fracInPeakGroup(:)))])
ylabel('Fraction of Correlograms with a given Peak Count')
legend(h,pkCountNames)
set(gca,'FontSize',plotProps.FigureFontSize)
title('Peak Counts Through Time')
box on
set(gca,'YTick',0:0.1:1)
ax = gca;
ax.YGrid = 'on';



%Plot the fraction of correlogram that are uniform and nonuniform
cc = pink(4); %set plot colors

figure()
%Indicate where injuries/treatments occurred
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2 2 -1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,max(fracInUnifGroup(:))+.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,max(fracInUnifGroup(:))+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
%Now plot the fraction of correlogram peaks in each classification
h=plot(1:AnalysisRegions.numRegions,fracInUnifGroup','.-','LineWidth',plotProps.lineWidth);
hold off
for c = 1:size(h,1)
    h(c).Color = cc(c,:);
    h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
end
set(gca,'XTick',1:AnalysisRegions.numRegions)
if isfield(AnalysisRegions,'RegionLabels')
    set(gca,'XTickLabel',AnalysisRegions.RegionLabels)
end

xlim([0.5 AnalysisRegions.numRegions+0.5])
xlabel('Recording Region')
ylim([0 1.1*(max(fracInUnifGroup(:)))])
ylabel('Fraction of Correlograms with a given Uniformity')
legend(h,unifNames)
set(gca,'FontSize',plotProps.FigureFontSize)
title('Correlogram Uniformity Through Time')
box on
set(gca,'YTick',0:0.1:1)
ax = gca;
ax.YGrid = 'on';

%Plot the fraction of correlogram peaks in each EEG classification
cc = pink(NumberOfFreqClassificationGroups+4); %set plot colors

figure()
%Indicate where injuries/treatments occurred
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[-1 2 2 -1],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,max(fracInFreqGroup(:))+.05,InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,max(fracInFreqGroup(:))+.05,[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
%Now plot the fraction of correlogram peaks in each classification
h=plot(1:AnalysisRegions.numRegions,fracInFreqGroup','.-','LineWidth',plotProps.lineWidth);
hold off
for c = 1:size(h,1)
    h(c).Color = cc(c,:);
    h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
end
set(gca,'XTick',1:AnalysisRegions.numRegions)
if isfield(AnalysisRegions,'RegionLabels')
    set(gca,'XTickLabel',AnalysisRegions.RegionLabels)
end
xlim([0.5 AnalysisRegions.numRegions+0.5])
xlabel('Recording Region')
ylim([0 1.2*(max(fracInFreqGroup(:)))])
if plotProps.classifyPeakTimesAsFrequencies == 1
    ylabel('Fraction of Correlogram Peaks in Each Frequency Band')
    title('Peak Frequency Classification Through Time')
else
     ylabel('Fraction of Correlogram Peaks in Each Time Range')
     title('Peak Location Classification Through Time')
end
legend(h,freqClassificationNames)
set(gca,'FontSize',plotProps.FigureFontSize)
box on
set(gca,'YTick',0:0.1:1)
ax = gca;
ax.YGrid = 'on';


%Plot the average plus/minus standard deviation of the frequency of each band
errorbarx = zeros(AnalysisRegions.numRegions,NumberOfFreqClassificationGroups); %set x array to give errorbar plots the correct x coordinate
for g = 1:NumberOfFreqClassificationGroups
    errorbarx(:,g)=1:AnalysisRegions.numRegions;
end

%%{
%error bars = range of the data
maxValy = maxFreqInFreqGroup-meanFreqInFreqGroup; maxVal = max(maxValy(:)); %find the maximum frequency
minValy = meanFreqInFreqGroup-minFreqInFreqGroup; minVal = min(min(minValy(:)),0.9); %find the minimum frequency
%}

%%{
%error bars = standard deviation of the data
maxValy = meanFreqInFreqGroup+sdFreqInFreqGroup; maxVal = max(maxValy(:)); %find the maximum frequency
minValy = meanFreqInFreqGroup-sdFreqInFreqGroup; minVal = min(min(minValy(:)),0.9); %find the minimum frequency
%}

%{
%error bars = standard error of the data
maxValy = meanFreqInFreqGroup+seFreqInFreqGroup; maxVal = max(maxValy(:)); %find the maximum frequency
minValy = meanFreqInFreqGroup-seFreqInFreqGroup; minVal = min(min(minValy(:)),0.9); %find the minimum frequency
%}

%get y-axis limits
if plotProps.classifyPeakTimesAsFrequencies == 1
    %Get powers of 10 to set y limits for the frequency plots
    vS = num2str(minVal);  b = strfind(vS,'0'); %find zeros in the minimum frequency
    if length(b)>1
        lowestPower = min(1+strfind(vS,'.'), b(2))-strfind(vS,'.'); %get the position after decimal where the minimum value starts (i.e. 1 = 0.1, 2 = 0.01, 3 = 0.001, 4 = 0.0001, etc...)
        minVal = 10.^-lowestPower;
    else %otherwise, mean - error went negative
        lowestPower = 1;
        minVal = 10.^-lowestPower;
    end

    highestPower = ceil(log10(maxVal)); %get the highest power to set upper y limit

else %otherwise classifying by peak times, not frequency
    lowestPower= -min(endpowr); 
    minVal = 0;
    highestPower = 0.1+max([stpowr endpowr]);
end

if plotProps.plotPeakLocationClassificationMedian == 1
    labelStr = 'Median';
else
    labelStr = 'Mean';
end

figure()
%Indicate where injuries/treatments occurred
hold on
if InjuryIndicies.NumberOfInjuriesOrTreatments > 0 %if there was an injury, shade where the injury occurred
    for j = 1:InjuryIndicies.NumberOfInjuriesOrTreatments %for all injuries/treatments...
        plot([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j)],[10.^(-lowestPower-1) 10.^highestPower],'--k','LineWidth',plotProps.lineWidth/2)
        plot([InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[10.^(-lowestPower-1) 10.^highestPower],'--k','LineWidth',plotProps.lineWidth/2)
        fill([InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)-InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j) InjuryInds(1,j)+InjuryInds(2,j)],[10.^(-lowestPower-1) 10.^highestPower 10.^highestPower 10.^(-lowestPower-1)],[0.6 0.6 0.6],'facealpha',0.2,'LineStyle','none')
        if InjuryInds(3,j) ~= 1 %if injury/treatment is mid-recoding
            text(InjuryInds(1,j) ,min(maxVal*1.5,9.5.^highestPower),InjuryIndicies.InjuryLabels{j},'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        else %otherwise injury/treatment is before recording start
            text(InjuryInds(1,j) ,min(maxVal*1.5,9.5.^highestPower),[InjuryIndicies.InjuryLabels{j} ' Before Recording Start'],'FontSize',plotProps.FigureFontSize*0.8,'HorizontalAlignment', 'center')
        end
    end
end
%Now plot mean plus/minus standard deviation in the frequencies of each group
%h =  errorbar(errorbarx,meanFreqInFreqGroup',minValy',maxValy','LineWidth',floor(plotProps.lineWidth*3/4));
if plotProps.plotPeakLocationClassificationMedian == 1
    h =  plot(errorbarx,medianFreqInFreqGroup','.-','LineWidth',floor(plotProps.lineWidth*3/4));
else
    h =  plot(errorbarx,meanFreqInFreqGroup','.-','LineWidth',floor(plotProps.lineWidth*3/4));
end
hold off
%for c = 1:size(h,2)
for c = 1:size(h,1)
    h(c).Color = cc(c,:);
    h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
end
set(gca,'XTick',1:AnalysisRegions.numRegions)
if isfield(AnalysisRegions,'RegionLabels')
    set(gca,'XTickLabel',AnalysisRegions.RegionLabels)
end
xlim([0.5 AnalysisRegions.numRegions+0.5])
xlabel('Recording Region')
ylim([minVal 10.^highestPower])
if plotProps.classifyPeakTimesAsFrequencies == 1
    ylabel([labelStr ' Frequency of Peaks in Each Band (Hz)'])
    title('Peak Frequencies Through Time')
else
    ylabel([labelStr ' Peak Times in Each Group (s)'])
    title('Peak Locations Through Time')
end
legend(h,freqClassificationNames)
set(gca,'FontSize',plotProps.FigureFontSize)
box on
set(gca,'YTick',10.^[-lowestPower:highestPower])
ax = gca;
ax.YGrid = 'on';
set(gca,'YScale','log')



%% Plot Classification Changes

%Plot group switching dynamics in different recording zones
cc = pink(numZones+2); %set plot colors

figure()
for g = 1:NumberOfleaderClassificationGroups %for each group...
    q = zeros(numZones,length(ZoneHistograms.LeaderFollower.Zone1));
    for zz = 1:numZones %cycle through zones...
        zoneName = strcat('Zone',num2str(zz)); %get the zone name to access the data structure
        histArray = ZoneHistograms.LeaderFollower.(zoneName); %for the zone, get the probabilities of being in current group, given past group

        q(zz,:) = ZoneHistograms.LeaderFollower.(zoneName)(g,:);

    end %end cycle through zones
    q = q';

    subplot(1,NumberOfleaderClassificationGroups,g)
    hold on
    h=plot(1:NumberOfleaderClassificationGroups, q,'-','MarkerSize',plotProps.markerSize,'LineWidth',plotProps.lineWidth);
    hold off
    for c = 1:size(h,1)
        h(c).Color = cc(c,:);
        h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
    end
    box on
    grid on
    title(['Previously Classified as ' leaderClassificationNames{g}])
    ylim([0 1])
    xlim([0 NumberOfleaderClassificationGroups+1])
    set(gca,'YTick',0:0.1:1)
    legend(groupOrder);
    ylabel('Probability of Current Classification')
    set(gca,'XTick',1:NumberOfleaderClassificationGroups,'XTickLabel',leaderClassificationNames)
    xlabel('Current Classification')
    set(gca,'FontSize',plotProps.FigureFontSize*0.7)

end


%Plot peak count switching dynamics in different recording zones
c = min(10,NumberOfpkCountClassificationGroups);
r = ceil(NumberOfpkCountClassificationGroups/c);


figure()
t=tiledlayout(r,c,'Padding','none','TileSpacing','compact','Padding','compact');
for g = 1:NumberOfpkCountClassificationGroups %for each group...
    q = zeros(numZones,length(ZoneHistograms.PeakCount.Zone1));
    for zz = 1:numZones %cycle through zones...
        zoneName = strcat('Zone',num2str(zz)); %get the zone name to access the data structure
        histArray = ZoneHistograms.PeakCount.(zoneName); %for the zone, get the probabilities of being in current group, given past group

        q(zz,:) = ZoneHistograms.PeakCount.(zoneName)(g,:);

    end %end cycle through zones
    q = q';

    nexttile
    hold on
    h=plot(1:NumberOfpkCountClassificationGroups, q,'-','LineWidth',plotProps.lineWidth);
    hold off
    for c = 1:size(h,1)
        h(c).Color = cc(c,:);
        h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
    end
    box on
    grid on
    title(['Previous Peak Count = ' pkCountNames{g}])
    ylim([0 1])
    xlim([0 NumberOfpkCountClassificationGroups+1])
    set(gca,'YTick',0:0.1:1)
    legend(groupOrder);
    set(gca,'XTick',1:NumberOfpkCountClassificationGroups,'XTickLabel',pkCountNames)
    set(gca,'FontSize',plotProps.FigureFontSize*0.7)

end
ylabel(t,'Probability of Current Peak Count','FontSize',plotProps.FigureFontSize)
xlabel(t,'Current Peak Count','FontSize',plotProps.FigureFontSize)


%Plot uniformity switching dynamics in different recording zones
c = min(10,NumberOfunifClassificationGroups);
r = ceil(NumberOfunifClassificationGroups/c);

figure()
t=tiledlayout(r,c,'Padding','none','TileSpacing','compact','Padding','compact');
for g = 1:NumberOfunifClassificationGroups %for each group...
    q = zeros(numZones,length(ZoneHistograms.Unif.Zone1));
    for zz = 1:numZones %cycle through zones...
        zoneName = strcat('Zone',num2str(zz)); %get the zone name to access the data structure
        histArray = ZoneHistograms.Unif.(zoneName); %for the zone, get the probabilities of being in current group, given past group

        q(zz,:) = ZoneHistograms.Unif.(zoneName)(g,:);

    end %end cycle through zones
    q = q';

    nexttile
    hold on
    h=plot(1:NumberOfunifClassificationGroups, q,'-','LineWidth',plotProps.lineWidth);
    hold off
    for c = 1:size(h,1)
        h(c).Color = cc(c,:);
        h(c).Marker = markerNmArr{c};   h(c).MarkerSize = markerSzArr(c);
    end
    box on
    grid on
    title(['Previously ' unifNames{g}])
    ylim([0 1])
    xlim([0 NumberOfunifClassificationGroups+1])
    set(gca,'YTick',0:0.1:1)
    legend(groupOrder);
    set(gca,'XTick',1:NumberOfunifClassificationGroups,'XTickLabel',unifNames)
    set(gca,'FontSize',plotProps.FigureFontSize*0.7)

end
ylabel(t,'Probability of Current Uniformity','FontSize',plotProps.FigureFontSize)
xlabel(t,'Current Uniformity','FontSize',plotProps.FigureFontSize)


%% Plot a summary of all possible correlogram classifications, and the percent of correlograms that fall into each category

totalCorrelogramCount = size(LFGroupAssignments,1)-sum(isnan(LFGroupAssignments)); %get total number of correlograms
UnifMat = zeros(NumberOfpkCountClassificationGroups,NumberOfleaderClassificationGroups,numZones); %initialize array to store classification combos of uniform correlograms
NonUnifMat = zeros(NumberOfpkCountClassificationGroups,NumberOfleaderClassificationGroups,numZones); %initialize array to store classification combos of nonuniform correlograms


for rr = 1:numZones %for each recording region or zone...
      rrInds = Zonez == rr;
    for u = 1:NumberOfunifClassificationGroups %for each uniformity classification...

        nm = strcat('G',num2str(u)); %get the uniformity group name
        UIdx=FindGroups.Unif.(nm)(:,rrInds); %find correlograms classified with the given value (nonuniform or uniform)

        for p = 1:NumberOfpkCountClassificationGroups %cycle through peak counts
            nm2 = strcat('G',num2str(p)); %get the peak count group name
            PIdx = FindGroups.PeakCount.(nm2)(:,rrInds); %find correlograms with the given peak count

            for l = 1:NumberOfleaderClassificationGroups %cycle through leader/follower groupings
                nm3 = strcat('G',num2str(l)); %get the leader/follower group name
                LIdx = FindGroups.LeaderFollower.(nm3)(:,rrInds); %find correlograms with the given leader/follower classification

                if u == 1 %check in which uniformity matrix to save the entry
                    NonUnifMat(p,l,rr) = sum(sum(UIdx&PIdx&LIdx))./sum(totalCorrelogramCount(rrInds));
                else %otherwise uniform matrix
                    UnifMat(p,l,rr) = sum(sum(UIdx&PIdx&LIdx))./sum(totalCorrelogramCount(rrInds));
                end %end save the fraction of correlograms with a given classification combo in the matrix

            end %end cycle through leader/follower groupings

        end %end cycle through peak counts

    end %end cycle through uniformity classifications

end %end cycle through zones

%Get x and y coordinates for the plot (these are peak counts and leader/follower classifications, uniformity is plotted as different markers)
x = pkCountGroupBoundaries;
y = 1:NumberOfleaderClassificationGroups;

colormapLength = 256;
cc = [1 1 1; jet(colormapLength)]; %get colormap to shade correlogram percentages
xposshift = 0;%0.1; %set a shift factor for the x position so uniform and uniform points are not plotted on top of eachother
yposshift = 0.07; %set a shift factor for the y position so uniform and uniform points are not plotted on top of eachother

for rr = 1:numZones %for each recording zone
    figure()
    hold on
    for xx = 1:length(x) %for each peak count...
        for yy = 1:length(y) %for each leader/follower class...
            %Plot uniform correlograms
            cmapIdx2 = max(1,1+ceil(colormapLength.*UnifMat(xx,yy,rr))); %shade the point according to % of correlograms with the given classification combo
            plot(x(xx)+xposshift,y(yy)+yposshift,'^','Color',cc(cmapIdx2,:),'MarkerSize',20,'MarkerFaceColor',cc(cmapIdx2,:),'MarkerEdgeColor',[0 0 0])
            text(x(xx)+xposshift+0.1,y(yy)+yposshift+0.1,[num2str(100.*UnifMat(xx,yy,rr),'%.1f') '%']) %display the % of correlograms with the given classification combo
            %Plot nonuniform correlograms
            cmapIdx = max(1,1+ceil(colormapLength.*NonUnifMat(xx,yy,rr))); %shade according to % of correlograms with the given classification combo
            plot(x(xx)-xposshift,y(yy)-yposshift,'v','MarkerSize',20,'MarkerFaceColor',cc(cmapIdx,:),'MarkerEdgeColor',[0 0 0])
            text(x(xx)-xposshift+0.1,y(yy)-yposshift-0.1,[num2str(100.*NonUnifMat(xx,yy,rr),'%.1f') '%']) %display the % of correlograms with the given classification combo
        end %end cycle through leader/follower classes
    end %end cycle through peak counts
     
    %Plot empty points for the legend
    h1=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength)),:),'MarkerEdgeColor',[0 0 0]);
    h2=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.9)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.9)),:),'MarkerEdgeColor',[0 0 0]);
    h3=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.8)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.8)),:),'MarkerEdgeColor',[0 0 0]);
    h4=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.7)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.7)),:),'MarkerEdgeColor',[0 0 0]);
    h5=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.6)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.6)),:),'MarkerEdgeColor',[0 0 0]);
    h6=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.5)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.5)),:),'MarkerEdgeColor',[0 0 0]);
    h7=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.4)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.4)),:),'MarkerEdgeColor',[0 0 0]);
    h8=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.3)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.3)),:),'MarkerEdgeColor',[0 0 0]);
    h9=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.2)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.2)),:),'MarkerEdgeColor',[0 0 0]);
    h10=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0.1)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0.1)),:),'MarkerEdgeColor',[0 0 0]);
    h11=plot(NaN,NaN,'s','Color',cc(max(1,round(colormapLength*0)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0)),:),'MarkerEdgeColor',[0 0 0]);

    hnu = plot(NaN,NaN,'v','Color',cc(max(1,round(colormapLength*0)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0)),:),'MarkerEdgeColor',[0 0 0]);
    hu = plot(NaN,NaN,'^','Color',cc(max(1,round(colormapLength*0)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0)),:),'MarkerEdgeColor',[0 0 0]);
    hblank = plot(NaN,NaN,'^','Color',cc(max(1,round(colormapLength*0)),:),'MarkerSize',20,'MarkerFaceColor',cc(max(1,round(colormapLength*0)),:),'MarkerEdgeColor',[1 1 1]);
    
    hold off
    %legend([hu,hnu,hblank, h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11], {'Uniform Correlograms';'Nonuniform Correlograms';'Color = Percent of all Correlograms:'; '100%';'90%';'80%';'70%';'60%';'50%';'40%';'30%';'20%';'10%';'0%'},'Location','EastOutside')
    legend([hu,hnu,h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11], {'Uniform Correlograms';'Nonuniform Correlograms';'100% Percent of all Correlograms';'90% Percent of all Correlograms';'80% Percent of all Correlograms';'70% Percent of all Correlograms';'60% Percent of all Correlograms';'50% Percent of all Correlograms';'40% Percent of all Correlograms';'30% Percent of all Correlograms';'20% Percent of all Correlograms';'10% Percent of all Correlograms';'0% Percent of all Correlograms'},'Location','EastOutside')
    box on
    grid on
    xlim([min(x)-10*max(xposshift,yposshift) max(x)+10*max(xposshift,yposshift)])
    ylim([min(y)-10*max(xposshift,yposshift) max(y)+10*max(xposshift,yposshift)])
    set(gca,'XTick',x,'XTickLabel',pkCountNames)
    xlabel('Number of Correlogram Peaks')
    set(gca,'YTick',y,'YTickLabel',leaderClassificationNames)
    ylabel('Leader/Follower Strength Classification')
    title(groupOrder{rr})
    set(gca,'FontSize',plotProps.FigureFontSize)
end %end cycle through recording zones


end