function Correlograms = generateCorrelograms(ProcessedData,InjuryIndicies,AnalysisRegions,sw,plotProps)
%This function makes correlograms by comparing each raster against all others
%       One raster acts as a reference, one is compared against the reference
%       All reference events are set as t=0 (one by one), the times of the comparison raster are shifted by t0, then a histogram is created using the resulting times
%       (note: only events that occurr within the correlogram bin width contribute to the correlogram for any given reference event)
%       Correlograms are analyzed for an entire analysis region (NOT the entire recording).  The number of correlograms made is the numberOfCells^2 (unless specific reference signals are selected, in which case it is numberOfCells x numberOfReferences) x numberOfRegions
%       Correlogram values are all normalized to probabilities (for ease of analysis later on), and the output sent to ProcessedData is the array of correlogram probabilities for each comparison

disp(' ')
disp('Generating correlograms for each of the selected regions.  This may take some time...')
disp(' ')

CorrelogramBins = linspace(-plotProps.correlogramBinMax,plotProps.correlogramBinMax,plotProps.numberOfCorrelogramBins); %make bins for the correlogram (these are bin CENTERS)
Correlograms.CorrelogramBins = CorrelogramBins; %save the bin centers in the correlogram structure

binSpacing = mean(gradient(CorrelogramBins)); %mean so not an entire array, just one number

BinEdges = sort([CorrelogramBins-(binSpacing/2) CorrelogramBins(end)+(binSpacing/2)]); %make correlogram bin edges


RegionTimes = AnalysisRegions.Seconds; %get the time cutoffs of each region
RegionTimes(RegionTimes==60) = 0; %if the region starts at the experiment start, set this to the first index (rather than 60 s, which was calculated from minute 1)

for rr = 1:AnalysisRegions.numRegions %cycle through regions of the recording being analyzed...

    %Get the portion of the rasters belonging to the given region
    InRegionInds = cellfun(@(x) find( RegionTimes(rr,1) <= x & x <= RegionTimes(rr,2) ),ProcessedData.Raster,'UniformOutput',false); %identify raster times that occurr within the given region (logical array)
    RegionRaster = cellfun(@(x,y) x(y),ProcessedData.Raster,InRegionInds,'UniformOutput',false); %Get the portion of the rasters that occurrs only in the given region (smaller raster containing only in-region events)

    NoFiringInTheRegionInds = cellfun(@isempty,RegionRaster); %check for signals without any events in the region and 
    RegionRaster(NoFiringInTheRegionInds) = {[NaN, NaN]'}; %Replace empty region rasters with NaN (so the empty values won't mess up the bsxfun portion of the code later, but also won't be part of a correlogram)

    correlogramIDs = cell(1,ProcessedData.CellCount^2); %initialize array to save signal names
    correlogramValues = zeros(length(CorrelogramBins),ProcessedData.CellCount^2); %initialize aorrelogram array
    numberOfSpikesContributingToTheCorrelogram = zeros(1,ProcessedData.CellCount^2); %initialize correlogram event count array

    columnIndex = 1; %initialize the variable directing where correlograms are saved
    for i = 1:ProcessedData.CellCount %cycle through reference cells...

        t0 = RegionRaster{i}; %get the reference times (all will be time 0 in the correlogram)

        RasterDiff = cellfun(@(x) bsxfun(@minus,x,t0'),RegionRaster,'UniformOutput',false); %calculate t - t0 for all rasters in the comparison and all t0s in the reference
        
        inCorrelogramInds = cellfun(@(x) x(x >= -plotProps.correlogramBinMax & x <= plotProps.correlogramBinMax) ,RasterDiff,'UniformOutput',false); %exclude t-t0 outside the correlogram bin range
        
        numSpikes = cellfun(@(x) length(x) ,inCorrelogramInds,'UniformOutput',false); %count the number of events from which the correlogram was constructed

        corrBinCounts = cellfun(@(x) histcounts(x,BinEdges,'Normalization','probability') ,inCorrelogramInds,'UniformOutput',false); %get the correlogram bin counts

        compNames = cellfun(@(x) [x  ' vs. ' ProcessedData.CellIDs{i}],ProcessedData.CellIDs,'UniformOutput',false); %get the name of the comparison vs. reference cell

        %correlogramTimez(:,columnIndex:columnIndex+ProcessedData.CellCount-1) = inCorrelogramInds; %save the times (that go into the correlogram)
        correlogramValues(:,columnIndex:columnIndex+ProcessedData.CellCount-1) = cell2mat(corrBinCounts')'; %convert the array containing the correlograms into a matrix
        correlogramIDs(columnIndex:columnIndex+ProcessedData.CellCount-1) = compNames; %enter the comparison list in the comparison array
        numberOfSpikesContributingToTheCorrelogram(1,columnIndex:columnIndex+ProcessedData.CellCount-1) = cell2mat(numSpikes); %enter the number of events that contributed to the correlogram into the array
        columnIndex = columnIndex+ProcessedData.CellCount; %update comparison column index for the next reference cell

    end %end cycle through reference cells

    %Save the correlogram of the given region
    reg_field = strcat('Region',num2str(rr)); %convert the name of the region to a variable name to be saved in the correlogram substructure
    Correlograms.(reg_field).CorrelogramProbabilities = correlogramValues; %save the correlogram (named by region) in the correlogram substructure
    Correlograms.(reg_field).ComparisonNames =  correlogramIDs; %save the names of the comparisons in each column of the correlogram
    Correlograms.(reg_field).numberOfSpikesContributingToTheCorrelogram = numberOfSpikesContributingToTheCorrelogram; %save the event count the contributed to each correlogram


    %Now make a matrix that saves indices where each signal is used in a correlogram
    %       row = the cell, column = indicies where the signal partakes in a correlogram as the comparison (CellAsComparisonInds matrix) or reference (CellAsReferenceInds matrix)

    comparisonStartIndex = strfind(correlogramIDs,'vs.'); %find where in the correlogram name "vs." begins (before = comparison cell)
    comparisonCellNames = cellfun(@(x,y) x(1:y-2) ,correlogramIDs, comparisonStartIndex ,'UniformOutput',false); %list only the name of the comparison cell
     
    refStartIndex = strfind(correlogramIDs,'.'); %find where in the correlogram name "vs." ends (after = reference cell)
    refCellNames = cellfun(@(x,y) x(y+2:end) ,correlogramIDs, refStartIndex ,'UniformOutput',false); %list only the name of the reference cell
     
    CellAsComparisonInds = zeros(ProcessedData.CellCount,ProcessedData.CellCount); %initialize indices matrix to save where the signal is used as a comparison
    CellAsReferenceInds = zeros(ProcessedData.CellCount,ProcessedData.CellCount); %initialize indices matrix to save where the signal is used as a reference

    for i = 1:ProcessedData.CellCount %cycle through cells...

        CellAsComparison = cell2mat(cellfun(@(x) contains(x,ProcessedData.CellIDs{i}), comparisonCellNames,'UniformOutput',false)); %find where the given signal was used as a comparison cell
        CellAsRef = cell2mat(cellfun(@(x) contains(x,ProcessedData.CellIDs{i}), refCellNames,'UniformOutput',false)); %find where the given signal was used as a reference cell

        CellAsComparisonInds(i,:) = find(CellAsComparison); %store the indices where the signal is used as a comparison in the comparison index matrix
        CellAsReferenceInds(i,:) = find(CellAsRef); %store the indices where the signal is used as a reference in the reference index matrix

    end %end cycle through cells

    %Save the index matricies
    %row = the cell, column = indicies where the signal partakes in a correlogram as the comparison (CellAsComparisonInds matrix) or reference (CellAsReferenceInds matrix)
    Correlograms.(reg_field).CellAsComparisonInds = CellAsComparisonInds;
    Correlograms.(reg_field).CellAsReferenceInds = CellAsReferenceInds;

end %end cycle through analysis regions


disp('     Done generating correlograms.')


end %end function