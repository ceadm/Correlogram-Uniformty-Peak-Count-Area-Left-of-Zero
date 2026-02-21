function plotCorrelograms(Correlograms,AnalysisRegions,RecordingMetrics,ProcessedData,plotProps)
%This function plots correlogram of specific comparison vs. reference signals that are specified by the user

donePlottingCorrelograms = 0; %initialize the logical indicating whether the user is finished plotting correlograms


%Label (or fetch labels) for each recording region
if ~isfield(AnalysisRegions,'RegionLabels') %if user automated region selection and hence did not type region labels...
    groupOrder = cell(1,AnalysisRegions.numRegions); %get each analysis region
    for rr = 1:AnalysisRegions.numRegions %for each region...
        groupOrder{rr} = 'Region' + " " + num2str(rr); %get the region name
    end %end cycle through regions

else %otherwise user assigned specific region labels, so use those instead
    groupOrder = AnalysisRegions.RegionLabels;
end


x = Correlograms.CorrelogramBins; %get the correlogram bin centers

%Get the number of rows and columns for the subplots of correlograms in each data analysis region
c = min(10,AnalysisRegions.numRegions); %number of columns for tile()
r = ceil(AnalysisRegions.numRegions/c); %number of rows for tile()

disp(' ')
disp(['Each recording region contains ' num2str(size(Correlograms.Region1.CorrelogramProbabilities,2)) ' correlograms.'])
disp('Type the name of the signals involved in the correlogram.  Note: these names should be obtained from heatmap points of interest.  Heatmaps can be displayed via switches in the main function.')
while donePlottingCorrelograms == 0
    CompLabel = input('     Please type the name of the COMPARISON signal. ','s');
    RefLabel = input('     Please type the name of the REFERENCE signal. ','s');

    %Make sure the names are lowercase so no capitalization issues ensue
    CompLabel = lower(CompLabel);
    RefLabel = lower(RefLabel);

    %Figure out the number assigned to the selected comparison signal
    compCellNumber = ismember(ProcessedData.CellIDs,CompLabel);
    refCellNumber = ismember(ProcessedData.CellIDs,RefLabel);

    %Find correlograms containing each signal in their respective roles
    idx1 = Correlograms.Region1.CellAsComparisonInds(compCellNumber,:);
    idx2 = Correlograms.Region1.CellAsReferenceInds(refCellNumber,:);

    if isempty(idx1) %the user mistyped the comparison name
        disp('The comparison signal entered was not part of the recording.  Please try again.')
    elseif isempty(idx2) %the user mistyped the reference name
        disp('The reference signal entered was not part of the recording.  Please try again.')

    else %As long as no mistakes were made in the signal names, get to plotting...

        idx = intersect(idx1,idx2); %find the particular correlogram

        %Get the max value of the correlogram to use in y-axis range
        M = 0; %max value in all correlograms of the given comparison vs. reference signals
        for rr = 1:AnalysisRegions.numRegions %for each region...

            reg_field = strcat('Region',num2str(rr)); %get the region name to access the data structur
            smoothedCorr = smoothdata(Correlograms.(reg_field).CorrelogramProbabilities(:,idx),'sgolay',plotProps.SmoothingWindowSize); %the smoothed correlogram used to calculate the number of peaks
            smoothedCorr = max(0,smoothedCorr); %make sure there are no negative entries
            smoothedCorr = smoothedCorr./sum(smoothedCorr); %normalize

            M = max(M,max(smoothedCorr)); %take the max of the max value of the smoothed correlogram used to calculate the number of peaks or a previous, higher value
        end


        %Get labels for the correlogram
        fullComp = char(Correlograms.(reg_field).ComparisonNames{idx});
        idx2 = cell2mat(strfind(Correlograms.(reg_field).ComparisonNames(idx),'vs.'));
        %Reassign comparison and reference names to make sure they match user input (reassignment only for debugging purposes)
        compName = fullComp(1:idx2-1);
        refName = fullComp(idx2+4:end);

        figure()

        t=tiledlayout(r,c,'Padding','none','TileSpacing','compact');
        for rr = 1:AnalysisRegions.numRegions %for each region...

            reg_field = strcat('Region',num2str(rr)); %get the region name to access the data structure
            smoothedCorr = smoothdata(Correlograms.(reg_field).CorrelogramProbabilities(:,idx),'loess',plotProps.SmoothingWindowSize);
            smoothedCorr = max(0,smoothedCorr); %make sure there are no negative entries
            smoothedCorr = smoothedCorr./sum(smoothedCorr); %normalize

            U=RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTest; %uniformity test decision
            Up=RecordingMetrics.CorrelogramMetrics.(reg_field).uniformityTestPValue; %uniformity test p value
            P=RecordingMetrics.CorrelogramMetrics.(reg_field).numberOfCorrelogramPeaks; %number of correlogram peaks
            L=RecordingMetrics.CorrelogramMetrics.(reg_field).leaderProb; %area left of zero

            nexttile
            plot([0 0],[-.9 M*1.9],'--k','LineWidth',plotProps.lineWidth/2)
            hold on
            bar(x,Correlograms.(reg_field).CorrelogramProbabilities(:,idx),'EdgeColor','None','FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(0.7,plotProps.faceAlpha))
            plot(x,smoothedCorr,'Color','k','LineWidth',plotProps.lineWidth)
            plot(RecordingMetrics.CorrelogramMetrics.(reg_field).correlogramPeakLocations{idx},RecordingMetrics.CorrelogramMetrics.(reg_field).correlogramPeaks{idx},'v','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'MarkerSize',plotProps.lineWidth*2)
            hold off
            text(0,M*1.6,['U = ' num2str(U(idx),'%.0f')],'HorizontalAlignment','center')
            text(0,M*1.4,['P = ' num2str(P(idx),'%.0f')],'HorizontalAlignment','center')
            text(0,M*1.2,['L = ' num2str(L(idx),'%.2f')],'HorizontalAlignment','center')
            ylim([0 M*1.8])
            xlim([min(x) max(x)])
            title([groupOrder{rr} num2str(Correlograms.(reg_field).numberOfSpikesContributingToTheCorrelogram(idx)) ' Events'])
            grid on
            box on
        end %end loop through regions
        ylabel(t,['Probability ' compName ' fired at time t given ' refName ' fired at t=0']);
        xlabel(t,'Time (s)');



        %Check whether the user wishes to plot more comparisons
        donePlottingCorrelograms =  input('Would you like to plot more correlograms?  Type 1 if yes, 0 if no. ');
        donePlottingCorrelograms = ~donePlottingCorrelograms;

    end %end check that both signals are part of the recording

end %end plotting correlograms

end %end function