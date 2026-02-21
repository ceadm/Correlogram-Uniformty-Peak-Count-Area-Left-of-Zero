function [RawData, AnalysisRegions, Correlograms, InjuryIndicies, ProcessedData, RecordingMetrics, plotProps] = analyzeRecordingPlxFile(plxFileName,sw,plotProps);
%This function initialized output folders and sequentially calls the
%functions to process and label the data, then generate and classify
%correlograms, then plot outputs

%% Initialize file paths

%Add folders to the current path to help matlab find the following...
%              The function to read .plx files is in 'Function to Read plx Files'
%              The functions that execute data analysis (for example, that generate rasters, correlograms, etc... are in 'Analysis Execution'
%              The folder that contains the raw data/.plx files is 'Raw Recordings'

%Check that all the necessary files are in place.  If not, the script will stop and tell what you need to obtain.
if ~exist('Function to Read plx Files','dir') || ~exist('Analysis Execution','dir') || ~exist('Raw Recordings','dir') || ~exist('Output Files','dir')

    if ~exist('Function to Read plx Files','dir'); disp('You need to obtain the read .plx function folder!  Approach a fellow labmate or go to https://www.mathworks.com/matlabcentral/fileexchange/42160-readplxfilec and click "Download" in the upper right.'); disp('     Function Citation: Benjamin Kraus (2024). readPLXFileC (https://www.mathworks.com/matlabcentral/fileexchange/42160-readplxfilec), MATLAB Central File Exchange.'); neededToInitialize = 1; end
    if ~exist('Analysis Execution','dir'); disp('You need to obtain the functions to excecute the code!  Obtain the folder from https://doi.org/10.1007/s12021-026-09770-9, or you will not be able to perform any calculations or analysis.'); neededToInitialize = 1; end
    if ~exist('Raw Recordings','dir'); mkdir('Raw Recordings\'); disp('Please put raw data (.plx files) into the "Raw Recordings" folder that was just created.'); neededToInitialize = 1; end
    if ~exist('Output Files','dir'); mkdir('Output Files\'); neededToInitialize = 0; end

else
    neededToInitialize = 0;

end

if neededToInitialize == 1
    disp('Code did not run because required folders are missing.  Please obtain the required items and try again.')

else  %if all the necessary folders are there, run the script

    %Make sure MATLAB knows that there are subfolders
    addpath([pwd '\Function to Read plx Files'],[pwd '\Analysis Execution'],[pwd '\Output Files']) %Add folders to the matlab path.  "pwd" is the current folder/directory

    %% Double check that parameters set by the user are within realistic ranges 
    
    %Ensure that plot settings are within acceptable ranges for the plot() function, and set to a default if not
    if plotProps.FigureFontSize < 1; plotProps.FigureFontSize = 10; disp('The desired figure font size was not a number greater than 0, so was reset to default.'); disp(' '); end %Plot font size should be larger than 0
    if plotProps.lineWidth < 1; plotProps.lineWidth = 2; disp('The desired figure line width was not a number greater than 0, so was reset to default.'); disp(' '); end  %Plot line width should be larger than 0
    if plotProps.markerSize < 1; plotProps.markerSize = 20; disp('The desired figure marker size was not a number greater than 0, so was reset to default.'); disp(' '); end %Plot marker size should be larger than 0
    if plotProps.faceAlpha < 0 || plotProps.faceAlpha > 1; plotProps.faceAlpha = 0.1; disp('The desired figure face alpha was not between 0 and 1, so was reset to default.'); disp(' '); end %Plot alpha should be between 0 and 1
    if size(plotProps.regionColors,2) ~= 3;  plotProps.regionColors = colorcube(10); disp('The desired figure color scheme did not have three columns, so was reset to default.'); disp(' '); end %Plot colors need to be an n x 3 array
    
    %Make sure that the automated region bin length is larger than 0, and set to a default if not
    if plotProps.automatedRegionBinLength <= 0; plotProps.automatedRegionBinLength = 1; disp('Automated region length was a number less than 0, so was set to 1 min.'); disp(' '); end 

    %Make sure that the p-value threshold to determine whether a correlogram is uniform or nonuniform is > 0, and set to 0.05 (default) if not
    if plotProps.unifPThresh <= 0; plotProps.unifPThresh = 0.05; disp('Uniformity cutoff p-value was set to 0 or negative, so was reset to default (0.05).'); disp(' '); end
    
    %Make sure that the correlogram time range is greater than 0, and set to a default of 1 s if not
    if plotProps.correlogramBinMax <= 0; plotProps.correlogramBinMax = 1; disp('Correlogram range was set to a time of 0 s or less, so was reset to default (1 s).'); disp(' '); end

    %Make sure that correlograms have at least 10 bins, and reset to a larger number (ensuring a bin size of 1 ms) if not
    if plotProps.numberOfCorrelogramBins < 10; plotProps.numberOfCorrelogramBins = ceil(plotProps.correlogramBinMax/.0005);  disp(['The desired number of correlogram bins was too small, so was reset to ' num2str(plotProps.numberOfCorrelogramBins) ' in order to create 1 ms bins throughout the specified correlogram range.']); disp(' '); end

    %Make sure that the max peak count for correlogram classifications is 2 or larger
    if plotProps.maxPeakCountBeforeNoise < 2; plotProps.maxPeakCountBeforeNoise = 10; disp('The desired max peak count used in correlogram classification was less than 2, so was reset to default (10).'); disp(' '); end %make sure that the max peak count before noise (for correlogram peak classification) is a number greater than 0


    %% Double check that the user has not set incompatible switches and that correlograms are centered around 0

    %Double check that the user has not entered an even bin number for correlograms (correlogram only has 0 centered if the bin number is odd, so odd bin counts must be ensured)
    if rem(plotProps.numberOfCorrelogramBins,2) == 0; plotProps.numberOfCorrelogramBins = plotProps.numberOfCorrelogramBins+1; end

     %If this is the first time running the script for a given recording, indicate that outputs (which don't yet exist) are not being loaded
    if ~exist([pwd '\Output Files\' plxFileName],'dir');  mkdir([pwd '\Output Files\' plxFileName]); sw.processDataFromScratch = 1; end

    %If the user wants to label bad signals or take only a segment of the recording, any pre-processing should be discarded, so the processDataFromScratch switch should be activated
    if sw.LabelJunkRegions == 1 || sw.LabelBadSignals == 1 || ~exist([pwd '\Output Files\' plxFileName '\badSignals.mat'],'file') || ~exist([pwd '\Output Files\' plxFileName '\badRegions.mat'],'file'); sw.processDataFromScratch = 1; sw.setAnalysisRegions = 1;  sw.setInjuryRegions = 1; end

    %If the user wants to set analysis regions, any pre-processing should be discarded, so the loadPreProcessedData switch should be inactivated
    if sw.setAnalysisRegions == 1; sw.processDataFromScratch = 1;  end

    %If this is the first time running the script, ensure that everything is calculated from scratch
    if ~exist([pwd '\Output Files\' plxFileName '\Correlograms.mat'],'file'); sw.processDataFromScratch = 1; end

    %If the stats are not saved, ensure that stats are calculated from scratch
    if ~exist([pwd '\Output Files\' plxFileName '\RecordingMetrics.mat'],'file'); sw.calculateStatisticsFromScratch = 1; end

    %If user changed correlogram properties, the correlograms and stats need to be regenerated.
    %Additionally, if recording region selection is different, correlograms need to be re-created and stats re-calculated
    if sw.processDataFromScratch == 0 && exist([pwd '\Output Files\' plxFileName '\plotProps.mat'],'file')
        OldPlotProps = load([pwd '\Output Files\' plxFileName '\plotProps.mat'],'plotProps');
        %Generate a logical to check whether the user changed any of the correlogram paramaters
        CorrelogramChangeIndicator = OldPlotProps.plotProps.automatedRegionBinLength ~= plotProps.automatedRegionBinLength || OldPlotProps.plotProps.unifPThresh ~= plotProps.unifPThresh || OldPlotProps.plotProps.correlogramBinMax ~= plotProps.correlogramBinMax || OldPlotProps.plotProps.numberOfCorrelogramBins ~= plotProps.numberOfCorrelogramBins || OldPlotProps.plotProps.correlogramSmoothingFactor ~= plotProps.correlogramSmoothingFactor || OldPlotProps.plotProps.sparseCorrelogramThresh ~= plotProps.sparseCorrelogramThresh;
        if CorrelogramChangeIndicator %if the user changed correlogram settings
            changedSettingsFlag = input('Correlograms are being loaded, but correlogram properties set in the main script differ from what was saved.  Would you like to remake the correlograms with the new properties?  Type 1 if yes, 0 otherwise.  If 0, the old correlograms will load. ');
            sw.processDataFromScratch = changedSettingsFlag; %the correlograms need to be regenerated and re-processed
           
            if changedSettingsFlag == 1 && OldPlotProps.plotProps.automatedRegionBinLength ~= plotProps.automatedRegionBinLength %if the region size has been altered...
                sw.setAnalysisRegions = 1; %re-assign recording regions
            end
            disp(' ')
        end
    end %end check whether user changed correlogram properties

    %Only load stats if the data is also being loaded.  Otherwise, stats need to be calculated fresh and this switch should be 1
    if sw.processDataFromScratch == 1;  sw.calculateStatisticsFromScratch = 1; end


    %% Set the desired figure font as default
    set(0, 'DefaultAxesFontName', plotProps.plotFont);
    set(0, 'DefaultUicontrolFontName', plotProps.plotFont);
    set(0, 'DefaultUitableFontName', plotProps.plotFont);
    set(0, 'DefaultAxesFontName', plotProps.plotFont);
    set(0, 'DefaultTextFontName', plotProps.plotFont);
    set(0, 'DefaultUipanelFontName', plotProps.plotFont);


    %% Load the Data
    RawData = readPLXFileC([pwd '\Raw Recordings\' plxFileName '.plx'],'fullread','spikes'); %load and read the .plx file
    %Function citation: Benjamin Kraus (2024). readPLXFileC (https://www.mathworks.com/matlabcentral/fileexchange/42160-readplxfilec), MATLAB Central File Exchange. Retrieved February 14, 2024.

    % ProcessedData is a structure containing processed data like rasters and correlograms.  It is initialized here in getRasters, where rasters are added to the structure
    ProcessedData = getRasters(RawData); % Generate Rasters (a.k.a. Spike Trains)
    % getRasters initializes ProcessedData and adds a list of signal IDs, a Raster (times when the given signal fired), and the total number of signals to ProcessedData

    %% Remove Signals that Died Out During the Recording or that were Noisy

    %If there are signals that the user wishes to exclude, remove those signals here
    ProcessedData = removeBadSignals(plxFileName,RawData,ProcessedData,sw);

    %If something happened during the recording, for example life support was being connected for the first m minutes of the recording, and you don't wish to include said data, you can eliminate that chunk of time in the next function
    ProcessedData = removeJunkRecordings(RawData,sw,plxFileName,ProcessedData);

    %% Plot and Calculate Raster Metrics
    %RecordingMetrics is a structure that contains any statistics on the data
    %       After calculateRasterMetrics, RecordingMetrics will contain the mean and standard deviation in events/min through time
    %calculateRasterMetrics takes the rasters in ProcessedData and makes:
    %       A plot of all rasters
    %       A histogram of the total number of times each signal fired
    %       A histogram of event intervals
    %       A plot of events/min vs. time
    [ProcessedData, RecordingMetrics] = calculateRasterMetrics(ProcessedData,plotProps);


    %% Locate Injuries or Treatments and Data Analysis (Correlogram) Regions
    %InjuryIndicies and AnalysisRegions are structures that give the time an injury or region to analyze begins (first column) and ends (second column). Units are in seconds or minutes as indicated in the name
    [InjuryIndicies,AnalysisRegions,plotProps] = getInjuryOrTreatmentIndicies(RecordingMetrics,plxFileName,sw,plotProps,ProcessedData);

    %Plot Regional Statistics
    ProcessedData = plotStatsInEachRegion(AnalysisRegions,InjuryIndicies,ProcessedData,plotProps,sw);

    %Pause to finish displaying the plots
    pause(1)

    %% Perform correlogram analysis

    if sw.processDataFromScratch == 1 || ~exist([pwd '\Output Files\' plxFileName '\Correlograms.mat'],'file') %if you are calculating correlograms from scratch and not loading previously-analyzed data...
        Correlograms = generateCorrelograms(ProcessedData,InjuryIndicies,AnalysisRegions,sw,plotProps); %create correlograms for each of the regions in AnalysisRegions
        save([pwd '\Output Files\' plxFileName '\Correlograms.mat'],'Correlograms')
    else %otherwise you want to load previously-generated correlograms
        load([pwd '\Output Files\' plxFileName '\Correlograms.mat'],'Correlograms');
    end

    if plotProps.correlogramSmoothingFactor <=0 %if the smoothing factor set by the user is negative or zero (won't work in the code), correct it to no smoothing
        plotProps.SmoothingWindowSize = 1; %update correlogram smoothing factor for peak counting
    else
        plotProps.SmoothingWindowSize = round(plotProps.correlogramSmoothingFactor*length(Correlograms.CorrelogramBins)); %update correlogram smoothing factor for peak counting
    end

    %% Get correlogram metrics and plots
    %Check whether metrics need to be recalculated or can just be loaded
    RecordingMetrics = calculateCorrelogramMetrics(Correlograms,AnalysisRegions,ProcessedData,RecordingMetrics,plotProps,sw, plxFileName);

    %Save all the matlab outputs
    if sw.processDataFromScratch == 1 %if re-processing or loading the data, save the new outputs
        save([pwd '\Output Files\' plxFileName '\RawData.mat'],'RawData')
        save([pwd '\Output Files\' plxFileName '\plotProps.mat'],'plotProps')
    end
    save([pwd '\Output Files\' plxFileName '\ProcessedData.mat'],'ProcessedData') %save this structure out of the statement because doing so doesn't add much time and keeps track of whether user wants general stats plotted by region or by treatment (otherwise changes to this setting wouldn't be saved)

    %If user wishes to plot specific correlograms from the data, do so
    if sw.plotIndividualCorrelograms == 1
        plotCorrelograms(Correlograms,AnalysisRegions,RecordingMetrics,ProcessedData,plotProps);
    end


    disp(' ')
    disp('Plotting correlogram statistics and classifications...')

    %Plot the outputs
    %Note all violin plots are generated using functions from https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot and Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847, https://www.mathworks.com/matlabcentral/fileexchange/170126-violinplot-matlab 
    PlotPopulationMetricChanges(ProcessedData,RecordingMetrics,AnalysisRegions,InjuryIndicies,Correlograms,plotProps,sw);

    %More Summary Plots
    summarizeCorrelograms(sw, plotProps, AnalysisRegions, InjuryIndicies, ProcessedData, RecordingMetrics, Correlograms);


end %end check that your folders are all there


end %end function