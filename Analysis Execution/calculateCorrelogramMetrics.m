function RecordingMetrics = calculateCorrelogramMetrics(Correlograms,AnalysisRegions,ProcessedData,RecordingMetrics,plotProps,sw, plxFileName)
%This function calculates correlogram metrics for each correlogram

if sw.calculateStatisticsFromScratch == 0 %if loading previously-calculated correlogram stats...
    disp(' ')
    disp('Loading correlogram statistics...')
    load([pwd '\Output Files\' plxFileName '\RecordingMetrics.mat'],'RecordingMetrics');

    %Identify specific variables in the structure (used later by the function for plotting)
    CorrelogramMetrics = RecordingMetrics.CorrelogramMetrics;
    RefLabels = CorrelogramMetrics.RefLabels;
    refIndicies = CorrelogramMetrics.refIndicies;
    maxPkCount = CorrelogramMetrics.maxPkCount;

else %otherwise calculate stats from scratch
    disp(' ')
    disp('Calculating correlogram statistics...')

    %If the user wishes to only use specific signals as the refrence signal
    if sw.SelectReferenceCells == 1 %if you wish to select specific reference signals (this feature allows the user to target specific subpopulations such as bursters)
        RefLabels = { }; %initialize reference list array
        refsDone = 0; %initialize switch to indicate you are still selecting reference signals
        disp(ProcessedData.CellIDs') %list signal IDs
        while refsDone == 0 %select reference signals...
            refLabel = input('Signal names are listed above.  Which would you like to use as the (next) reference signal? Type the signal name here.  ','s');
            RefLabels = [RefLabels, refLabel];
            refsDone = input('     Are you finished selecting reference signals?  Type 1 if yes, 0 if no. ');
        end %end select reference signals

        RefLabels = unique(RefLabels); %make sure there are no duplicate names in the reference list
        refIndicies = zeros(1,length(RefLabels)); %initialize a list to find reference signal indicies
        for r = 1:length(RefLabels) %for each selected reference...
            b=contains(ProcessedData.CellIDs,RefLabels(r)); %find the signal ID that matches the reference ID
            refIndicies(r) = find(b); %get the index of the reference signal
        end %end get reference index
    else
        refIndicies = 1:ProcessedData.CellCount; %otherwise use every signal as a reference
        RefLabels = ProcessedData.CellIDs; %get all signal labels for figures
    end

    CorrelogramMetrics.RefLabels = RefLabels; %save the reference signal IDs in the correlogram metric structure
    CorrelogramMetrics.refIndicies = refIndicies; %save the reference signal indicies in the correlogram metric structure


    %The following metrics are calculated for individual correlograms

    % Uniformity quantifies whether there is a relation between the signals (if the distribution is uniform, firing of the two signals is independent)
    %      Determined using a chi^2 test (applies to binned distributions
    %      regardless of how the binning is performed)

    % The number of correlogram peaks quantifies the firing signature of the signals (more peaks = more complicated pattern)

    % Correlogram area left or right of zero quantifies the leader or follower character of the signal
    %      If more area is left of zero, the comparison signal has a tendency to proceed firing of the reference
    %      If more area is right of zero, the comparison signal has a tendency to follow the reference
    %      If the area right and left are similar, there is no consistent order

    leftOfZeroInds = Correlograms.CorrelogramBins<0; %find correlogram bins that occurr before t=0
    rightOfZeroInds = Correlograms.CorrelogramBins>0; %find correlogram bins that occurr after t=0
    zeroInds = Correlograms.CorrelogramBins==0; %find the correlogram bin at t=0
    rightOfZeroInds(end) = []; %the last is a bin edge that should be removed so correlograms and rightOfZeroInds have the same length
    
    maxPkCount = 0; %initialize variable to store the maximum number of correlogram peaks for all correlograms (used later for plotting)

    for rr = 1:AnalysisRegions.numRegions %for each analyzed region with correlograms...

        reg_field = strcat('Region',num2str(rr)); %convert the name of the region to a variable name to call the correct correlogram substructure
        rcor = Correlograms.(reg_field).CorrelogramProbabilities; %get the correlograms of the given region
        compInds = Correlograms.(reg_field).CellAsComparisonInds; %get indicies of where each signal is used as a comparison
        refInds = Correlograms.(reg_field).CellAsReferenceInds; %get indicies of where each signal is used as a reference

        %find correlograms with too few events/sparse distributions based on user threshold
        PotentiallySparse = Correlograms.(reg_field).numberOfSpikesContributingToTheCorrelogram < plotProps.sparseCorrelogramThresh;

        %Get leader/follower character
        leaderProb = sum(rcor(leftOfZeroInds,:),1) + sum(rcor(zeroInds,:),1)./2; %get area left of zero
        followerProb = sum(rcor(rightOfZeroInds,:),1) + sum(rcor(zeroInds,:),1)./2; %get area right of zero (note: leaderProb + followerProb = 1)

        %Generate a uniform array to compare against correlograms to perform a statistical test for uniformity
        unifarray = ones(size(Correlograms.CorrelogramBins)); unifarray = unifarray./sum(unifarray);

        %First, smooth the correlograms in order to limit detection of false peaks or false uniformity
        %smoothedCorr = smoothdata(rcor,'gaussian',plotProps.SmoothingWindowSize); %works, but can sometimes smooth out refractory periods
        smoothedCorr = smoothdata(rcor,'loess',plotProps.SmoothingWindowSize); %less likely to smooth out refractory periods
        smoothedCorr = max(0,smoothedCorr); %make sure there are no negative entries (some smoothing algorithms allow for negatives, so ensure that doesn't happen here)
        smoothedCorr = smoothedCorr./sum(smoothedCorr,1); %normalize to probability

        %Test for uniform distributions
        %Since correlograms are binned/categorical data (not continuous), a chi squared test is used to test for uniformity.  If p < user-specified significance cutoff, then the correlogram is not uniform
        N = Correlograms.(reg_field).numberOfSpikesContributingToTheCorrelogram; %get the number of events in the correlogram to calculate the chi^2 value - N is used to convert from probability to count to avoid smaller numbers and thereby improve test accuracy
        unifTestp = cell2mat(arrayfun(@(col)  chi2cdf(sum((((N(col).*rcor(:,col))'-(N(col).*unifarray)).^2)./(N(col).*unifarray)),length(Correlograms.CorrelogramBins)-1,'upper'), 1:size(rcor,2),'UniformOutput',false)); %get the probability that the distribution is uniform using a chi^2 test
        unifTest = unifTestp > plotProps.unifPThresh; %the chi2 test is not significant if uniform and correlogram distributions are the same.  Therefore, to mark uniformity, find where the chi2 test output was not significant

        %Get the peak values (pks) and the time bin indicies at which the peaks occurr (loc) for each correlogram
        %Also, ensures peak prominence should be above the uniform array prominence so tiny fluctuations are not counted as peaks
        [pks,loci] = arrayfun(@(col)findpeaks(smoothedCorr(:,col),'MinPeakProminence',unifarray(1)),1:size(smoothedCorr,2),'UniformOutput',false);

        %Convert peak indicies to peak times (s)
        loc = cellfun(@(x) Correlograms.CorrelogramBins(x),loci,'UniformOutput',false);

        %Count the number of correlogram peaks
        NumPeaks = cellfun(@length,pks); %get the number of peaks in each correlogram

        %Set sparse correlograms to NaN because one cannot trust that these values are not due to having too few events in the raster
        leaderProb(PotentiallySparse) = NaN;
        followerProb(PotentiallySparse) = NaN;
        %unifTest(PotentiallySparse) = NaN; %Here for display only.  NaN values cannot be converted to logicals, so unifTest = 0 can either be because correlogram is empty/has too few vents or because the correlogram is nonuniform
        pks(PotentiallySparse) = {[NaN]};
        loc(PotentiallySparse) = {[NaN]};
        NumPeaks(PotentiallySparse) = NaN;

        %Save all the metrics for the given region in the correlogram metrics structure
        CorrelogramMetrics.(reg_field).leaderProb = leaderProb;
        CorrelogramMetrics.(reg_field).followerProb = followerProb;
        CorrelogramMetrics.(reg_field).uniformityTest = unifTest;
        CorrelogramMetrics.(reg_field).uniformityTestPValue = unifTestp;
        CorrelogramMetrics.(reg_field).correlogramPeaks = pks;
        CorrelogramMetrics.(reg_field).correlogramPeakLocations = loc;
        CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks = NumPeaks;


        %Now initialize matricies for plotting correlogram metrics (rows = comparison signal, columns = reference signal)
        leaderProbMat = zeros(ProcessedData.CellCount,ProcessedData.CellCount); %area left of zero matrix
        unifTestMat = zeros(ProcessedData.CellCount,ProcessedData.CellCount); %uniformity matrix
        unifTestpMat = zeros(ProcessedData.CellCount,ProcessedData.CellCount); %p matrix (from chi^2 test, used to determine unifTest decision)
        numPeaksMat = zeros(ProcessedData.CellCount,ProcessedData.CellCount); %number of correlogram peaks matrix
        
        for i = 1:ProcessedData.CellCount %cycle through signals so each signal as a comparison gets the relevant matrix row...
            leaderProbMat(i,:) = leaderProb(compInds(i,:)); %make a matrix to visualize probability of being a leader
            unifTestMat(i,:) = unifTest(compInds(i,:)); %make a matrix to visualize uniformity
            unifTestpMat(i,:) = unifTestp(compInds(i,:)); %make a matrix to visualize uniformity probability
            numPeaksMat(i,:) = NumPeaks(compInds(i,:)); %make a matrix to visualize the number of correlogram peaks
        end %end cycle through signals

        %Save the matricies created in the previous loop
        CorrelogramMetrics.(reg_field).leaderProbMat = leaderProbMat;
        CorrelogramMetrics.(reg_field).numPeaksMat = numPeaksMat;
        CorrelogramMetrics.(reg_field).unifTestMat = unifTestMat;
        CorrelogramMetrics.(reg_field).unifTestpMat = unifTestpMat;
        CorrelogramMetrics.maxPkCount = maxPkCount;

        maxPkCount = max([maxPkCount max(numPeaksMat(:))]); %update the maximum peak count (of all correlograms)

    end %end cycle through analysis regions

   
end %end check whether to load the metrics from a previous run or calculate metrics from scratch

%% Now plot the correlogram metrics and calculate differences in the metrics between regions

for rr = 1:AnalysisRegions.numRegions %for each analyzed region with correlograms...

    reg_field = strcat('Region',num2str(rr)); %convert the name of the region to a variable name to call the correct correlogram substructure

    if sw.plotCorrelogramLeaderProbabilities == 1 %if you want to plot leader probabilities
        %Plot leadership probability (area left of 0)
        Q = CorrelogramMetrics.(reg_field).leaderProbMat(:,refIndicies); Q(isnan(Q)) = -1000; %convert NaNs to low number for plot coloration
        figure()
        imagesc(Q)
        axis square
        box on
        set(gca,'XTick',1:length(refIndicies),'XTickLabel',RefLabels)
        set(gca,'YTick',1:ProcessedData.CellCount,'YTickLabel',ProcessedData.CellIDs)
        caxis([-.01 1])
        colormap([1 1 1; jet(256)])
        cbr = colorbar;
        tl = get(cbr, 'TickLabels');
        tl{1} = 'Not Enough Events in the Correlogram'; % blank range indicator as string
        set(cbr, 'TickLabels', tl)
        ylabel('Comparison Cell')
        xlabel('Reference Cell')
        set(gca,'FontSize',plotProps.FigureFontSize/2)
        if ~isfield(AnalysisRegions,'RegionLabels')
            title(['Correlogram Area Left of Zero in Recording Region ' num2str(rr)])
        else
            title(['Correlogram Area Left of Zero ' AnalysisRegions.RegionLabels{rr}])
        end

    end %end check to plot leader probabilities

    %Plot independence (whether distribution is uniform or not)
    if sw.plotCorrelogramUniformities == 1 %if you want to plot correlogram uniformity
        Q = CorrelogramMetrics.(reg_field).unifTestMat(:,refIndicies);
        Q(isnan(Q)) = -1000; %convert NaNs to low number for plot coloration
        figure()
        imagesc(Q) %uniformityMat
        axis square
        box on
        caxis([-1 1])
        set(gca,'XTick',1:length(refIndicies),'XTickLabel',RefLabels)
        set(gca,'YTick',1:ProcessedData.CellCount,'YTickLabel',ProcessedData.CellIDs)
        colormap([1 1 1; 0 .2 0; 0 1 0])
        cbr = colorbar('YTick',[-1 0 1],'YTickLabel', num2cell([-1 0 1]));
        tl = get(cbr, 'TickLabels');
        tl{1} = 'Not Enough Events in the Correlogram'; % blank range indicator as string
        tl{2} = 'Nonuniform Distribution';
        tl{3} = 'Uniform Distribution';
        set(cbr, 'TickLabels', tl)
        ylabel('Comparison Cell')
        xlabel('Reference Cell')
        set(gca,'FontSize',plotProps.FigureFontSize/2)
        if ~isfield(AnalysisRegions,'RegionLabels')
            title(['Uniform Correlograms in Recording Region ' num2str(rr)])
        else
            title(['Uniform Correlograms ' AnalysisRegions.RegionLabels{rr}])
        end

        %Now plot chi^2 p used in uniformity test decision
        Q = CorrelogramMetrics.(reg_field).unifTestpMat(:,refIndicies);
        Q = log10(Q); %convert to log10 scale for ease of viewing
        Q(Q<log10(eps)) = log10(eps); %if Q lower than machine error (eps), set to machine error
        Q(isnan(Q)) = -1000; %convert NaNs to low number for plot coloration
        figure()
        imagesc(Q) %uniformitypMat
        axis square
        box on
        caxis([round(log10(eps))-.01 0])
        colormap([1 1 1; jet(256)])
        cbr = colorbar;
        tl = get(cbr, 'TickLabels');
        tl{1} = 'Not Enough Events in the Correlogram'; % blank range indicator as string
        set(cbr, 'TickLabels', tl)
        ylabel('Comparison Cell')
        xlabel('Reference Cell')
        set(gca,'FontSize',plotProps.FigureFontSize/2)
        if ~isfield(AnalysisRegions,'RegionLabels')
            title(['log_{10}(Probability of a Correlogram being Uniform in Recording Region ' num2str(rr) ')'])
        else
            title(['log_{10}(Probability of a Correlogram being Uniform ' AnalysisRegions.RegionLabels{rr} ')'])
        end

    end %end check whether to plot correlogram uniformity

    %Plot number of correlogram peaks
    if sw.plotCorrelogramNumberOfPeaks == 1 %if you want to plot peak counts
        Q = CorrelogramMetrics.(reg_field).numPeaksMat(:,refIndicies); Q(isnan(Q)) = -1000;
        figure()
        imagesc(Q)
        axis square
        box on
        caxis([0 maxPkCount])
        set(gca,'XTick',1:length(refIndicies),'XTickLabel',RefLabels)
        set(gca,'YTick',1:ProcessedData.CellCount,'YTickLabel',ProcessedData.CellIDs)
        colormap([1 1 1; copper(256)])
        cbr = colorbar;
        tl = get(cbr, 'TickLabels');
        tl{1} = 'Not Enough Events in the Correlogram'; % blank range indicator as string
        set(cbr, 'TickLabels', tl)
        ylabel('Comparison Cell')
        xlabel('Reference Cell')
        set(gca,'FontSize',plotProps.FigureFontSize/2)
        if ~isfield(AnalysisRegions,'RegionLabels')
            title(['Number of Correlogram Peaks in Recording Region ' num2str(rr)])
        else
            title(['Number of Correlogram Peaks ' AnalysisRegions.RegionLabels{rr}])
        end

    end %end check whether to plot peak counts



    %Now compare the metrics between regions and plot the difference (plot only if desired)
    if rr>1 %if you are not on the first region...
        for rr2 = 1:rr-1 %cycle through regions before the current region (after will be compared to the current rr the next time rr advances)...

            reg_field2 = strcat('Region',num2str(rr2)); %convert the name of the region to a variable name to call the correct correlogram substructure

            diffname = strcat(reg_field,'_Minus_'); diffname = strcat(diffname,reg_field2); %get a variable name for the difference between the two regions

            %Get difference in leader probability
            Q = CorrelogramMetrics.(reg_field).leaderProbMat(:,refIndicies) - CorrelogramMetrics.(reg_field2).leaderProbMat(:,refIndicies); %area left of zero in current - previous region
            CorrelogramMetrics.RegionDifferences.(diffname).leaderProbDifference = Q; %save the leader probability difference matrix
            Q(isnan(Q)) = -1000; %convert NaNs to low number for plot coloration
            if sw.plotDifferencesInCorrelogramMetricsBetweenRegions == 1 && sw.plotCorrelogramLeaderProbabilities == 1 %if you want to plot the differences between regions...
                figure()
                imagesc(Q)
                axis square
                box on
                set(gca,'XTick',1:length(refIndicies),'XTickLabel',RefLabels)
                set(gca,'YTick',1:ProcessedData.CellCount,'YTickLabel',ProcessedData.CellIDs)
                caxis([-1 1])
                colormap([1 1 1; jet(256)])
                colorbar('YTick',[-1:0.1:1],'YTickLabel', [{'Not Enough Events in the Correlogram'}, num2signal(-.9:0.1:1)])
                ylabel('Comparison Cell')
                xlabel('Reference Cell')
                set(gca,'FontSize',plotProps.FigureFontSize/2)
                title({'Change in Correlogram Area Left of Zero'; ['(Region ' num2str(rr) ' - Region ' num2str(rr2) ')']})
                if ~isfield(AnalysisRegions,'RegionLabels')
                    title({'Change in Correlogram Area Left of Zero'; ['(Region ' num2str(rr) ' - Region ' num2str(rr2) ')']})
                else
                    title({'Change in Correlogram Area Left of Zero'; ['(' AnalysisRegions.RegionLabels{rr} ' - ' AnalysisRegions.RegionLabels{rr2} ')']})
                end
            end %end check if want to plot differences in the metrics between each region

            %Plot independence changes (whether distribution is uniform or not)
            Q = CorrelogramMetrics.(reg_field).unifTestMat(:,refIndicies) - CorrelogramMetrics.(reg_field2).unifTestMat(:,refIndicies); %uniformity of present - previous regions
            CorrelogramMetrics.RegionDifferences.(diffname).uniformityDifference = Q; %save the uniformitydifference matrix
            Q(isnan(Q)) = -1000; %convert NaNs to low number for plot coloration
            if sw.plotDifferencesInCorrelogramMetricsBetweenRegions == 1  && sw.plotCorrelogramUniformities == 1 %if you want to plot the differences between regions...
                figure()
                imagesc(Q) %uniformity difference
                axis square
                box on
                caxis([-2 1])
                set(gca,'XTick',1:length(refIndicies),'XTickLabel',RefLabels)
                set(gca,'YTick',1:ProcessedData.CellCount,'YTickLabel',ProcessedData.CellIDs)
                colormap([1 1 1; hsv(3)])
                cbr = colorbar('YTick',[-2 -1 0 1],'YTickLabel', num2signal([-2 -1 0 1]));
                set(cbr, 'TickLabels', {'Not Enough Events in the Correlogram','Lost Uniformity','No Uniformity Change','Gained Uniformity'})
                ylabel('Comparison Cell')
                xlabel('Reference Cell')
                set(gca,'FontSize',plotProps.FigureFontSize/2)
                if ~isfield(AnalysisRegions,'RegionLabels')
                    title({'Change in Correlogram Uniformity'; ['(Region ' num2str(rr) ' - Region ' num2str(rr2) ')']})
                else
                    title({'Change in Correlogram Uniformity'; ['(' AnalysisRegions.RegionLabels{rr} ' - ' AnalysisRegions.RegionLabels{rr2} ')']})
                end

            end %end check if want to plot differences in the metrics between each region


            %Plot firing signature (number of correlogram peaks)
            Q = CorrelogramMetrics.(reg_field).numPeaksMat(:,refIndicies) - CorrelogramMetrics.(reg_field2).numPeaksMat(:,refIndicies); %number of correlogram peaks in the current - previous region
            CorrelogramMetrics.RegionDifferences.(diffname).peakCountDifference = Q; %save the peak count difference matrix
            lowestQ = min(Q(:)); highestQ = max(Q(:)); Qlims = max(abs([lowestQ highestQ])); %get the min and max for plot scaling later
            Q(isnan(Q)) = -1000; %set NaN to low value for ease of coloring
            if sw.plotDifferencesInCorrelogramMetricsBetweenRegions == 1 && sw.plotCorrelogramNumberOfPeaks == 1 %if you want to plot the differences between regions...
                figure()
                imagesc(Q)
                axis square
                box on
                caxis([-Qlims  Qlims])
                set(gca,'XTick',1:length(refIndicies),'XTickLabel',RefLabels)
                set(gca,'YTick',1:ProcessedData.CellCount,'YTickLabel',ProcessedData.CellIDs)
                cc1 = copper(256/2); cc2 = copper(256/2); cc2(:,2)=0; cc = [flipud(cc1); (cc2)]; colormap([1 1 1; cc])
                cbr = colorbar;
                tl = get(cbr, 'TickLabels');
                tl{1} = 'Not Enough Events in the Correlogram'; % blank range indicator as string
                set(cbr, 'TickLabels', tl)
                ylabel('Comparison Cell')
                xlabel('Reference Cell')
                set(gca,'FontSize',plotProps.FigureFontSize/2)
                if ~isfield(AnalysisRegions,'RegionLabels')
                    title({'Change in the Number of Correlogram Peaks'; ['(Region ' num2str(rr) ' - Region ' num2str(rr2) ')']})
                else
                    title({'Change in the Number of Correlogram Peaks'; ['(' AnalysisRegions.RegionLabels{rr} ' - ' AnalysisRegions.RegionLabels{rr2} ')']})
                end

            end %end check if want to plot differences in the metrics between each region
        end %end cycle through comparison regions
    end %end check that the loop is not on the first region

end %end cycle through analysis regions


%Save the Output
RecordingMetrics.CorrelogramMetrics = CorrelogramMetrics;
save([pwd '\Output Files\' plxFileName '\RecordingMetrics.mat'],'RecordingMetrics')


end %end function
