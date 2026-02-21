function ProcessedData = removeBadSignals(plxFileName,a,ProcessedData,sw)
%This function allows the user to exclude specific signals from analysis if
%desired.  It will run through a series of prompts for the user to list the
%signal(s) to remove.

if sw.LabelBadSignals == 1 || ~exist([pwd '\Output Files\' plxFileName '\badSignals.mat'],'file') %if you want to eliminate signals or this is the first time running the script for the given recording...

    %% Remove signals that died during the recording or contained too much noise
    badSignalsIndicator = input('Are there any signals you wish to exclude?  Type 1 if yes, 0 otherwise. ');

    if badSignalsIndicator == 1

        NoiseLabels = { }; %initialize list of signals to ignore
        noiseDone = 0; %initialize switch to indicate you are still selecting signals to remove
        disp(' ')
        disp('List of Active Signals in the Recording: ')
        for i = 1:ProcessedData.CellCount
            disp(['    ' ProcessedData.CellIDs{i} ' (during the experiment, this signal was labeled ' ProcessedData.CellDisplayNames{i} ')'])
        end

        while noiseDone == 0 %select reference signals...
            NoiseLabel = input(' Signal names are listed above.  Which are noisy or lost during recording? Type the name here.  ','s');
            if ~any(ismember(ProcessedData.CellIDs,NoiseLabel))  && ~any(ismember(ProcessedData.CellDisplayNames,NoiseLabel)) %if the user typed something that is not a signal name
                disp('     Error:  The name entered is not a signal in this recording.  Please try again. ')

            else %otherwise, the user entered a legit name
                NoiseLabels = [NoiseLabels, NoiseLabel];
                noiseDone = input('     Are you finished listing signals to exclude?  Type 1 if yes, 0 if no. ');
            end

        end %end select noisy or lost signals

        %Make sure none of the noisy signals are listed twice
        NoiseLabels = unique(NoiseLabels);
        ProcessedData.CellIDs
        ProcessedData.CellDisplayNames
        %Now create a logical array indicating which signals to exclude
        badSignalInds = cellfun(@(a) find(ismember(ProcessedData.CellIDs,a) | ismember(ProcessedData.CellDisplayNames,a)),NoiseLabels,'UniformOutput',false); %find signals with data that is not empty

        badSignals = false(1,ProcessedData.CellCount); %initialize logial array to indicate signals to ignore
        for j = 1:length(badSignalInds) %cycle through all bad signals...
            badSignals(badSignalInds{j})=true; %mark which signals are bad in the logical array
        end %end cycle through bad signals

    else %otherwise all signals are good
        badSignals = false(1,ProcessedData.CellCount);

    end %end check whethere there are any signals to exclude

    save([pwd '\Output Files\' plxFileName '\badSignals.mat'],'badSignals')


else %otherwise load pre-assigned signal exclusion indicies

    load([pwd '\Output Files\' plxFileName '\badSignals.mat'],'badSignals');

end %end check whether to exclude signals (set or loaded)

%Now remove the desired signals from the data
ProcessedData.CellIDs(badSignals) =[]; %remove from signal name list
ProcessedData.Raster(badSignals) =[]; %remove from raster list
ProcessedData.CellCount = length(ProcessedData.CellIDs); %get the updated signal count

Unitz = {a.SpikeChannels.Units}; %tells which unit(s) fired at each timestamp in TS
%Find electrodes that had one or more signals
activeSignals = ~cellfun(@isempty,Unitz); %find units with data that is not empty

ProcessedData.NumberOfActiveElectrodes = sum(activeSignals); %save the number of active electrodes
ProcessedData.TotalNumberOfElectrodes = length(activeSignals); %save the total number of electrodes

end %end function