function RasterStructure = getRasters(a)
%Reads a plx file and obtains the recording raster and signals

%Save any comments in the file
RasterStructure.Comment = a.Comment;

%Find the total gain
RasterStructure.ExperimentTotalGain = [a.SpikeChannels.Gain].*a.SpikePreAmpGain;

TS = {a.SpikeChannels.Timestamps}; %tells the time at which firing events occurred for each signal/cell
TS = cellfun(@(x) double(x)./a.ADFrequency,TS,'UniformOutput',false); %convert from timestamp to seconds

Unitz = {a.SpikeChannels.Units}; %tells which cells/signals on a given electrode (4 cells max) fired at each time in TS

%Mark locations that had one or more active signals
activeSignals = ~cellfun(@isempty,Unitz); %find electrodes that are not empty

%Find the electrode to which each signal belongs
electrodeID = cell2mat({a.SpikeChannels.SIG});
displayName = {a.SpikeChannels.SIGName};

cellCount = 0; %initialize signal counter

for i = find(activeSignals) %cycle through electrodes...

    q = unique(Unitz{i}); %find the number of signals on that electrode

    spikeTimez = TS{i}; %get the times when any signal fired on the electrode

    for j = 1:length(q) %cycle through all signals on that particular electrode...

        idx = Unitz{i}==q(j); %find events belonging to the given signal only

        cellRaster = spikeTimez(idx); %get the times of events in the raster

        %Now get the signal name
        dspNm = displayName{i}; %Get the name of the signal as displayed during the recording (note: this name is independent of electrode on which the signal is located)
        numberStart = strfind(dspNm,'0'); numberStart = max(numberStart(numberStart~=length(dspNm))); %remove zeros padding the numbers (this turns 002 into 2).  However, if there is a zero at the end of the number (i.e. 010, 020, 030...), ensure that zero does not count (i.e. so 020 correctly becomes 20). 
        if j == 1 %if on the first cell
            cellName = ['s' num2str(electrodeID(i)) 'a']; %signal name consists of three parts: (1) s for "signal", followed by (2) a number indicating the electrode at which the signal is located, followed by (3) the letter "a" (indicating this is the first signal on the electrode)
            cellDisplayName = 'dsp' + convertCharsToStrings(dspNm(1+numberStart:end)) + 'a'; %signal name as displayed during the recording (for reference purposes)
        elseif j == 2 %if on the second cell
            cellName = ['s' num2str(electrodeID(i)) 'b'];  %signal name consists of three parts: (1) s for "signal", followed by (2) a number indicating the electrode at which the signal is located, followed by (3) the letter "b" (indicating this is the second signal on the electrode)
            cellDisplayName = 'dsp' + convertCharsToStrings(dspNm(1+numberStart:end)) + 'b'; %signal name as displayed during the recording (for reference purposes)
        elseif j == 3 %if on the third cell
            cellName = ['s' num2str(electrodeID(i)) 'c'];  %signal name consists of three parts: (1) s for "signal", followed by (2) a number indicating the electrode at which the signal is located, followed by (3) the letter "c" (indicating this is the third signal on the electrode)
            cellDisplayName = 'dsp' + convertCharsToStrings(dspNm(1+numberStart:end)) + 'c'; %signal name as displayed during the recording (for reference purposes)
        elseif j == 4 %if on the fourth cell
            cellName = ['s' num2str(electrodeID(i)) 'd'];  %signal name consists of three parts: (1) s for "signal", followed by (2) a number indicating the electrode at which the signal is located, followed by (3) the letter "d" (indicating this is the fourth signal on the electrode)
            cellDisplayName = 'dsp' + convertCharsToStrings(dspNm(1+numberStart:end)) + 'd'; %signal name as displayed during the recording (for reference purposes)
        else
            disp(['More than four signals on electrode ' num2str(electrodeID(i)) '.  Any more than four per electrode cannot be analyzed.'])
        end

        cellCount = cellCount+1; %update the signal count

        %Now store the signal name and raster
        RasterStructure.CellIDs{cellCount} = cellName; %save the signal name in a name array
        RasterStructure.CellDisplayNames{cellCount} = char(cellDisplayName); %save the signal name in a name array
        RasterStructure.Raster{cellCount} = cellRaster;  %save the raster in a raster array

    end %end cycle through cells

end %end cycle through electrodes

RasterStructure.CellCount = cellCount;


end %end function


