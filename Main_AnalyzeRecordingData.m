set(0,'DefaultFigureWindowStyle','Docked')
%close all hidden; close all; 
clear; 
clc

tic %start a timer to see how long the algorithm runs

%This script loads a .plx file, makes rasters from the plx, makes correlograms from the rasters, then analyzes data from the rasters and correlograms
%Key metrics are plotted.  If you wish to save the plots, you must do so manually
%All the .m files (matlab files generated in the course of the analysis) are saved in the 'Outputs' folder

%NOTE ON SAVING OUTPUTS LARGER THAN 2GB
%If a warning occurrs that a file was too large to save, the file with name FileName can be saved manually via the command window by typing: save([pwd '\Output Files\' plxFileName '\FileName.mat'],'FileName','-v7.3')
%    For example, depending on recording length and region size, the variable "Correlograms" can be large, and MATLAB's default saving will send a warning.  In this event, Correlograms can manually be saved by typing save([pwd '\Output Files\' plxFileName '\Correlograms.mat'],'Correlograms','-v7.3') in the command window

%% User Controlled Parameters
plxFileName = 'ER52_ImpactWithBicuculline';
% The plx files used in the original publication are located in "Raw Recordings", and are named as follows.
%     'SMJM_Bicuculline' = recording of a network treated with 40 uM bicuculline methiodide
%     'SM_pHshock' = recording of a network subject to CO2 shutoff and subsequent alkalosis
%     'ER52_ImpactWithBicuculline' = recording of a network under bicuculline subjected to impact injury
%                This recording was previously published, and reused with permission 
%                Original Publication: Rogers, E. A., & Gross, G. W. (2019). Simultaneous electrophysiological and morphological assessment of functional damage to neural networks in vitro after 30–300 g impacts. Scientific Reports, 9(1), 14994. https://doi.org/10.1038/s41598-019-51541-x
%     'JMSM_ImpactWithoutBicuculline' = recording of a network NOT under bicuculline subjected to impact injury and subsequent control treatments (details in the paper's SI)

%% Switches turn optional features of the algorithm on and off.  The switch is on when set to 1
%%%%%%%% Switches for data reprocessing (ignored if this is the first time the given recording is being analyzed) %%%%%%%% 
sw.LabelJunkRegions = 0; %set to 1 if you wish to indicate whether and where there are regions of the recording to discard (for example, due to a power outage, equipment failure, etc...), 0 loads pre-processed data
sw.LabelBadSignals = 0; %set to 1 if you wish to label noisy or lost signals, 0 if you want to load previously-indicated bad signals
sw.setAnalysisRegions = 0; %set to 1 if you wish to select recording regions and re-analyze the data, 0 if you wish to load previously-selected regions
sw.setInjuryRegions = 0; %set to 1 if you wish to reselect injury/treatment times, 0 to load previously-set injury/treatment regions from a previous analysis

sw.processDataFromScratch = 0; %set to 1 if you want to process data (make correlograms) from scratch, 0 to load previously analyzed data (if you have previously processed the data and just want to skip to plotting)
sw.calculateStatisticsFromScratch = 0; %set to 1 if you want to calculate correlogram stats from scratch, 0 to load previously calculated stats (if you have previously processed the data and just want to skip to plotting)

%%%%%%%% These switches are helpful for limiting the number of plots (of correlogram metrics) generated at once %%%%%%%% 
sw.plotCorrelogramNumberOfPeaks = 0; %set to 1 if you want to plot correlogam peak count for each region of the recording being analyzed
sw.plotCorrelogramUniformities = 0; %set to 1 if you want to plot correlogam uniformity for each region of the recording being analyzed
sw.plotCorrelogramLeaderProbabilities = 0; %set to 1 if you want to plot correlogam area left of zero (leader probability) for each region of the recording being analyzed
sw.plotDifferencesInCorrelogramMetricsBetweenRegions = 0; %set to 1 if you want to plot differences in correlogram metrics between each region of the recording 

sw.plotIndividualCorrelograms = 0; %set to 1 if you wish to plot the correlograms of specific signal pairings

sw.plotRecordingStatsForEachDataRegion = 0; %set to 1 if you wish to plot event counts and event intervals for each recording region 
%Note: if above switch is set to 0 and there is an injury or treatment during the recording, then the script will make plots for before vs. after treatment.  For example, if there are two treatments, there will be three histograms for: (1) before anything, (2) after treatment 1, but before treatment 2, (3) after treatment 2

plotProps.numberOfRegionsForASinglePlot = 1e9; %set the number of regions that can be plotted on the same graph.  If more regions than this parameter, matlab will use a tiled figure layout with one tile = one region.  If there are fewer or equal regions to this parameter, they will all be plotted on a single plot

%% plotProps is a structure that contains properties for plot visualization and correlogram properties

%%%%%%%% Figure appearance %%%%%%%% 
plotProps.FigureFontSize = 18; %sets the font size for all plots 
plotProps.lineWidth = 2; %sets the line width for all plots 
plotProps.markerSize = 20; %sets the marker size for all plots 
plotProps.faceAlpha = 0.1; %sets shading transparency for all plots 
plotProps.regionColors = colorcube(10); %sets violin plot and region colors
plotProps.plotFont = 'Arial'; %sets the font for all plots 

%%%%%%%% Correlogram properties %%%%%%%% 
plotProps.automatedRegionBinLength = 7; %units = minutes, if you do not wish to manually select regions of the data to analyze, the code will automatically divide the recording into x minute intervals.  x is set with this parameter here

plotProps.unifPThresh = 0.05; %p-value below which correlograms are considered nonuniform
plotProps.correlogramBinMax = 1; %time around t0 you want to examine in seconds (for example 0.2 here will mean the correlogram's x axis will range from -0.2 s to 0.2 s)
plotProps.numberOfCorrelogramBins = ceil(plotProps.correlogramBinMax/.0005); %this is the number of bins in the correlogram - if you prefer to specify by bin width, divide binMax by half the desired bin width and round up
plotProps.correlogramSmoothingFactor = 0.01; %for correlogram peak counting: this is a factor you multiply by correlogram bin length to smooth the data to reduce false peak detections (smaller = less smoothing)
plotProps.sparseCorrelogramThresh = 0; %this is the minimum number of events that contribute to a correlogram for the correlogram to be considered trustworthy.  This switch allows the user to ignore sparse correlograms to avoid misinterpretating the data.
%Note: if the above threshold is zero, there may be some correlograms with 0 events (if a signal never fired in the given region).  Such correlograms will show as white dots in the matrix plots controlled by sw

plotProps.maxPeakCountBeforeNoise = 10; %threshold above which the peak count is grouped into one classification
plotProps.classifyPeakTimesAsFrequencies = 0; %set to 1 if you want to classify peak times (t) by frequency (f = 1/t, Hz).  If 0, peak times will be classified as times (s).
plotProps.plotPeakLocationClassificationMedian = 0; %set to 1 if you want to plot the median peak location value (f or t) within each classification bin.  If set to 0, the mean will be plotted instead of the median.

plotProps.includeAutocorrelogramsInStatistics = 0; %set to 1 if autocorrelograms should be included in recording statistics (otherwise, only crosscorrelograms will be considered)

sw.SelectReferenceCells = 0; %set to 1 if you would like to hand-select reference signals and use only a subset of the popoulation as references.  If set to 0, all signals are used as references

%% Perform data analysis with the parameters specified above
addpath([pwd '\Analysis Execution']) %Add folder to the matlab path so matlab can find analysis functions.  "pwd" is the current folder/directory

[RawData, AnalysisRegions, Correlograms, InjuryIndicies, ProcessedData, RecordingMetrics, plotProps] = analyzeRecordingPlxFile(plxFileName,sw,plotProps); %analye data of the given recording

disp(' ') %display a space before displaying elapsed time
toc %end the timer and display script runtime