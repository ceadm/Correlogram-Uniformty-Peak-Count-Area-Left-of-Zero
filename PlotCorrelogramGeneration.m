%This script was used to make figure 1 of 

%% TO RUN THIS SCRIPT...
%  FIRST RUN Main_AnalyzeRecordingData.m ON SMJM_Bicuculline.  
%  THEN RUN THIS SCRIPT WITHOUT CLEARING ANYTHING.

%%

%Below is copy/pasted from generateCorrelograms.m, except for the plotting.

CorrelogramBins = linspace(-plotProps.correlogramBinMax,plotProps.correlogramBinMax,plotProps.numberOfCorrelogramBins); %make bins for the correlogram (these are bin CENTERS)
Correlograms.CorrelogramBins = CorrelogramBins; %save the bin centers in the correlogram structure

binSpacing = mean(gradient(CorrelogramBins)); %mean so not an entire array, just one number

BinEdges = sort([CorrelogramBins-(binSpacing/2) CorrelogramBins(end)+(binSpacing/2)]);


RegionTimes = AnalysisRegions.Seconds; %get the time cutoffs of each region
RegionTimes(RegionTimes==60) = 0; %if the region starts at the experiment start, set this to the first index (rather than 60 s, which was calculated from minute 1)


correlogramValues = zeros(length(CorrelogramBins),1);

for rr = 9%1:AnalysisRegions.numRegions %cycle through regions of the recording being analyzed...

    %Get the portion of the rasters belonging to the given region
    InRegionInds = cellfun(@(x) find( RegionTimes(rr,1) <= x & x <= RegionTimes(rr,2) ),ProcessedData.Raster,'UniformOutput',false); %identify raster times that occurr within the given region
    RegionRaster = cellfun(@(x,y) x(y),ProcessedData.Raster,InRegionInds,'UniformOutput',false); %Get the portion of the rasters that occurrs only in the given region

    CompRaster = RegionRaster{strcmp(ProcessedData.CellIDs,'s39a')}; %comparison cell

    for i = find(strcmp(ProcessedData.CellIDs,'s30a'))%1:ProcessedData.CellCount %cycle through reference cells...

        t0 = RegionRaster{i}; %get the reference times (all will be time 0 in the correlogram)

        for rs = 1:length(t0) %for each ref spike...

            ct = CompRaster - t0(rs); %comparison spike time

            ct = ct(abs(ct)<= max(BinEdges));

            corrBinCounts = histcounts(ct,BinEdges); %get the correlogram bin counts

            correlogramValues = corrBinCounts' + correlogramValues;

            if ismember(rs,[1 2 10 40 200 round(length(t0)/2) length(t0)])  %for particular reference spikes (tailored to the specific region, so adjust as needed)

                smoothedCorr = smoothdata(correlogramValues,'sgolay',plotProps.SmoothingWindowSize); %the smoothed correlogram used to calculate the number of peaks
                smoothedCorr = max(0,smoothedCorr); %make sure there are no negative entries
                %smoothedCorr = smoothedCorr./sum(smoothedCorr); %normalize

                startTime = find(ct>=min(CorrelogramBins),1,'first');
                endTime = find(ct<=max(CorrelogramBins),1,'last');

                startTimeR = find(t0-t0(rs)>=min(CorrelogramBins),1,'first');
                endTimeR = find(t0-t0(rs)<=max(CorrelogramBins),1,'last');

                figure() %plot the correlogram at the current reference spike

                subplot(1,3,1) %plot the current section of the raster
                plot(ct(startTime:endTime)+t0(rs),zeros(size(ct(startTime:endTime)))+2,'|k','LineWidth',plotProps.lineWidth/2,'MarkerSize',plotProps.markerSize*4)
                hold on
                plot(t0(startTimeR:endTimeR),zeros(size(t0(startTimeR:endTimeR)))+1,'|','Color',[0.7 0.7 0.7],'LineWidth',plotProps.lineWidth/2,'MarkerSize',plotProps.markerSize*4)
                plot(t0(rs),1,'|m','LineWidth',plotProps.lineWidth,'MarkerSize',plotProps.markerSize*4)
                hold off
                legend('Comparison','Reference', ['Reference Spike ' num2str(rs)])
                xlabel('Spike Times (s)')
                ylabel('Raster')
                box on; grid on;
                xlim([t0(rs)-1 t0(rs)+1])
                set(gca,'YTick',[1 2],'YTickLabel',{'Reference','Comparison'})
                ylim([0 3])
                title({['Raster Events Near Reference Spike ' num2str(rs)];['(' num2str(sum(corrBinCounts>0)) ' Comparison Spikes)']})
                set(gca,'FontSize',plotProps.FigureFontSize)

                subplot(1,3,2) %plot the current contrubution to the correlogram
                bar(CorrelogramBins,corrBinCounts','EdgeColor',0.6.*plotProps.regionColors(rr,:),'FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(0.7,plotProps.faceAlpha))
                title({'Contribution to Correlogram';['(' num2str(sum(corrBinCounts>0)) ' Entries)']})
                xlabel('Time (s)')
                ylabel('Count')
                box on; grid on
                set(gca,'FontSize',plotProps.FigureFontSize)
                set(gca,'XTick',-1:.1:1,'XTickLabel',{'-1','','-0.8','','-0.6','','-0.4','','-0.2','','0','','0.2','','0.4','','0.6','','0.8','','1'}); xtickangle(0)

                subplot(1,3,3) %plot the overall correlogram
                hold on
                bar(CorrelogramBins,correlogramValues,'EdgeColor','None','FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(1,plotProps.faceAlpha))
                plot(CorrelogramBins,smoothedCorr,'Color','k','LineWidth',plotProps.lineWidth)
                hold off
                legend('Correlogram','Smoothed Correlogram')
                title({['Correlogram After Reference Spike ' num2str(rs)]; ['(' num2str(sum(correlogramValues)) ' Entries Total)']})
                xlabel('Time (s)')
                ylabel('Count')
                set(gca,'XTick',-1:.1:1,'XTickLabel',{'-1','','-0.8','','-0.6','','-0.4','','-0.2','','0','','0.2','','0.4','','0.6','','0.8','','1'}); xtickangle(0)
                ylim([0 max(correlogramValues)+1])
                box on; grid on
                set(gca,'FontSize',plotProps.FigureFontSize)
            end

        end %end cycle through ref spikes

        figure() %plot the final correlogram
        subplot(1,3,3)
        hold on
        bar(CorrelogramBins,correlogramValues./sum(correlogramValues),'EdgeColor','None','FaceColor',plotProps.regionColors(rr,:),'FaceAlpha',max(1,plotProps.faceAlpha))
        plot(CorrelogramBins,smoothedCorr./sum(smoothedCorr),'Color','k','LineWidth',plotProps.lineWidth)
        hold off
        legend('Correlogram','Smoothed Correlogram')
        title('Normalized Correlogram')
        xlabel('Time (s)')
        ylabel('Probability')
        set(gca,'XTick',-1:.1:1,'XTickLabel',{'-1','','-0.8','','-0.6','','-0.4','','-0.2','','0','','0.2','','0.4','','0.6','','0.8','','1'}); xtickangle(0)
        ylim([0 max(correlogramValues./sum(correlogramValues))])
        box on; grid on
        set(gca,'FontSize',plotProps.FigureFontSize)

    end %end cycle through ref



end %end cycle through analysis regions


