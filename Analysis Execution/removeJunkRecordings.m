function ProcessedData = removeJunkRecordings(a,sw,plxFileName,ProcessedData)
% This function allows the user to delete recording segments that should not be part of the anaysis
% For example, if life support was being connected during the first few minutes of the recording,
% or if life support failed at the recording end,
% and the user doesn't want to include the affected recording segments

%Save the true experiment end and recording sampling frequency
ProcessedData.ExperimentEndTimestamp = a.LastTimestamp; %save the experiment end in the array
ProcessedData.ExperimentSamplFreq = a.ADFrequency; %save the experiment sampling frequency so you can convert the timestamps to seconds, minutes, hours, days, whatever you want

%Get the experiment end time in seconds (the recording ran from second 0 to this end time)
ExperimentEndTime = a.LastTimestamp/a.ADFrequency; %time the experiment ended in seconds

if sw.LabelJunkRegions == 1 || ~exist([pwd '\Output Files\' plxFileName '\badRegions.mat'],'file') %if you want to enter bad regions of the data or this is the first time running the script for the given recording...
    %% Label bad regions
    %Logical indicating whether the user had the recordings described above
    badRegions.badRegionsIndicator = input('Are there sections of the recording you wish to exclude from analysis?  Type 1 if yes, 0 otherwise.  Note: such sections should have been noted during the recording or must be examined in the recording software before running this script. ');

    if badRegions.badRegionsIndicator == 1 %if there are regions that recorded data with confounding variables...

        RegionsEliminated = 0; %logical indicating whether the user is done typing the start and end times of regions to eliminate
        badRegions.badStartTime = []; badRegions.badEndTime = []; %initialize arrays to save the start and end of recording sections to delete
        i = 1; %initialize counter for number of segments to delete
        while RegionsEliminated == 0 %while labeling bad regions...
            badRegions.badStartTime(i) = input('     What minute did the section to delete start? If experiment start, type 0. ');
            badRegions.badEndTime(i) = input(['     What minute did the section to delete end? If experiment end, type ' num2str(ceil(ExperimentEndTime/60)) '. ']);

            badRegionMistakeIndicator = input('     Was there a mistake in any of the entries for this section?  Type 1 if yes, 0 if no.  If yes, do not fret.  You will re-enter the correct information for this section. '); 
            if badRegionMistakeIndicator ~= 1 %if no mistakes were made...
                i = i+1; %advance bad region counter
            end
            RegionsEliminated = input('Are you finished listing recording sections to delete?  Type 1 if yes, 0 otherwise. ');
        end %end label bad regions

        %Raster times are all in seconds.  Therefore, convert bad start/end time lists from minutes to seconds so units match
        badRegions.badStartTime = badRegions.badStartTime.*60;
        badRegions.badEndTime = badRegions.badEndTime.*60;

    end %end check whether to label junk signals

    save([pwd '\Output Files\' plxFileName '\badRegions.mat'],'badRegions')

else %otherwise load pre-labeled junk recording regions

    load([pwd '\Output Files\' plxFileName '\badRegions.mat'],'badRegions');

end %end check whether to load junk indicies

%% Delete Bad Regions
if badRegions.badRegionsIndicator == 1 %if there are regions of the data to remove...
    %Delete bad regions
    for i = 1:length(badRegions.badStartTime) %cycle through bad regions...
        if i == 1 %for the first region to delete, initialize a logical array marking events in each raster that fall in the regions to delete
            spikesInBadRegion  = cellfun(@(x) (x>=badRegions.badStartTime(i) & x<=badRegions.badEndTime(i)),ProcessedData.Raster,'UniformOutput',false);
        else  %add to the logical array marking events in each raster that fall in the regions to delete
            spikesInBadRegion = cellfun(@(x,y) y | (x>=badRegions.badStartTime(i) & x<=badRegions.badEndTime(i)),ProcessedData.Raster,spikesInBadRegion,'UniformOutput',false);
        end
    end %end cycle through bad regions

    %Now remove any events in these regions from consideration in the raster
    b = cellfun(@(x,y) x(~y),ProcessedData.Raster,spikesInBadRegion,'UniformOutput',false); %list only the events that occurred outside of weird regions in the recording

    %If several minutes of the recording end were eliminated, update the recording end time
    if any(badRegions.badEndTime>ExperimentEndTime) %if part of the recording end is being eliminated...
        %Set the actual recording end to the start time of the trimmed end region minus 1 ms or the max value left in the rasters
        ExperimentEndTime = max(badRegions.badStartTime(badRegions.badEndTime>ExperimentEndTime)-0.001,  max(max(cell2mat(cellfun(@(x) max(x),b,'UniformOutput',false)'))+0.001));
    end %end check whether recording start is being eliminated

    %If several minutes of the recording start were eliminated, shift the times (so that everything starts at time 0 again)
    if any(badRegions.badStartTime==0) %if part of the recording start is being eliminated...
        %Set the new experiment start to either: 1 ms before the first event OR the end time of the cut region (whichever is smaller)
        minz = min(min(cell2mat(cellfun(@(x) min(x),b,'UniformOutput',false)'))-0.001, badRegions.badEndTime(badRegions.badStartTime==0));

        %Subtract this start time from the event times
        b = cellfun(@(x) x-minz,b,'UniformOutput',false); %list only the events that occurred outside of weird regions in the recording
        ExperimentEndTime = ExperimentEndTime-minz;
    end %end check whether recording start is being eliminated

    %Determine whether any of the rasters are now empty
    b2 = cell2mat(cellfun(@(x) isempty(x),b,'UniformOutput',false));

    ProcessedData.Raster = b(~b2);
    ProcessedData.CellIDs = ProcessedData.CellIDs(~b2);
    ProcessedData.CellCount = sum(~b2);

end %end check for regions of the data to remove

ProcessedData.ExperimentEndTime = ExperimentEndTime; %save the experiment end time (in s) in the array

end