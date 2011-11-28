function UTIL_plotMatrix(toPlot, method)
% UTIL_plotMatrix general purpose plotting utility for plotting the results
%  of the MAP auditory model.
% All plots are placed in subplots of a figure (default figure 1).
%
% Input arguments:
% 	'toPlot' is matrix (either numeric or logical)
% 	'method' is a structure containing plot instructions
%
% mandatory parameters:
% 	method.displaydt		xValues spacing between data points
%   method.yValues          yaxis labels mandatory only for 3D plots
%
% optional
% 	method.figureNo         default figure(1)
% 	method.numPlots         number of subPlots in the figure (default=1)
% 	method.subPlotNo        number of this plot (default=1)
% 	method.zValuesRange     [min max] value pair to define yaxis limits
%   method.zValuesRange     [min max] CLIMS for 3-D plot
% 	method.yLabel           (string) y-axis label
%   method.minyMaxy         y-axis limits
% 	method.xLabel           (string) x-axis label
% 	method.title  		    (string) subplot title
%   method.bar    		    =1,  to force bar histogram (single channel only)
%   method.view             3D plot 'view' settings e.g. [-6 40]
%   method.axes             (handle) where to plot (overules all others)
%   method.maxPixels        maximum number of pixels (used to speed plotting)
%   method.blackOnWhite     =1; inverts display for 2D plots
%   method.forceLog         positive values are put on log z-scale
%   method.rasterDotSize    min value is 1
%   method.defaultFontSize  deafult= 12
%   method.timeStart        default= dt
%   method.defaultTextColor default ='k'
%   method.defaultAxesColor default ='k'
%   method.nCols            default = 1 (layout for subplots)
%   method.nRows            default=method.numPlots (layout for subplots
%   method.segmentNumber    plot only this segment while 'hold on'
%
% e.g.
%   UTIL_plotMatrix(toPlot, method)

dt=method.displaydt;
[r cols]=size(toPlot);
if cols==1
    % toPlot should be a wide matrix or a long vector
    toPlot=toPlot';
end

if ~isfield(method,'numPlots') || isempty(method.numPlots)
    method.numPlots =1;
    method.subPlotNo =1;
end

if ~isfield(method,'figureNo') || isempty(method.figureNo)
    method.figureNo=99;
end

% if ~isfield(method,'zValuesRange') || isempty(method.zValuesRange)
%     method.zValuesRange=[-inf inf];
% end

% set some defaults
if ~isfield( method,'blackOnWhite') || isempty(method.blackOnWhite)
    method.blackOnWhite=0;
end
if ~isfield(method,'timeStart')|| isempty(method.timeStart)
    method.timeStart=dt;
end
if ~isfield(method,'objectDuration') || isempty(method.objectDuration)
    [nRows nCols]=size(toPlot); method.objectDuration=dt*nCols;
end
if ~isfield(method,'defaultFontSize') || isempty(method.defaultFontSize)
    method.defaultFontSize=12;
end
if ~isfield(method,'defaultTextColor') || isempty(method.defaultTextColor)
    method.defaultTextColor='k';
    defaultTextColor=method.defaultTextColor;
else
    defaultTextColor='k';
end
if ~isfield( method,'defaultAxesColor') || isempty(method.defaultAxesColor)
    method.defaultAxesColor=defaultTextColor;
end
defaultAxesColor=method.defaultAxesColor;

% arrangement of plots in rows and columns
if ~isfield(method,'nCols') || isempty(method.nRows)
    method.nCols=1;
end
if  ~isfield(method,'nRows') || isempty(method.nRows)
    method.nRows= method.numPlots;
end

if ~isfield(method,'rasterDotSize') || isempty(method.rasterDotSize)
    rasterDotSize=1;
else
    rasterDotSize=method.rasterDotSize;
end

% user can specify either an independent axis
%   or a subplot of the current figure
%   if both are specified, 'axes' takes priority
figure(method.figureNo)
if isfield(method,'axes') && ~isempty(method.axes)
    % user defines where to plot it
    axes(method.axes);
    method.numPlots =1;
    method.subPlotNo =1;
    
else
    % now using a regular figure
    if method.subPlotNo>method.numPlots;
        error('UTIL_plotMatrix: not enough subplots allocated in figure 1.  Check method.numPlots')
    end
    % choose subplot
    subplot(method.nRows,method.nCols,method.subPlotNo),  % cla
    
    if isfield(method,'segmentNumber') && ~isempty(method.segmentNumber)...
            && method.segmentNumber>1
        % in multi-segment mode do not clear the image
        %  from the previous segment
        hold on
    else
        % otherwise a fresh image will be plotted
        hold off
        cla
    end
end

[numYvalues numXvalues]=size(toPlot);
xValues=method.timeStart:dt:method.timeStart+dt*(numXvalues-1);

if isfield(method,'yValues') && ~isempty(method.yValues)
    % yValues is normally a vector specifying channel BF
    yValues=method.yValues;
else
    yValues=1:numYvalues;
end

if round(numYvalues/length(yValues))>1
    % case where the plot matrix is double height (e.g. LSR+HSR)
    yValues=[yValues yValues];
    method.plotDivider=1;
else
    method.plotDivider=0;
end

% Now start the plot.
%  3D plotting for 4 or more channels
%  otherwise special cases for fewer channels

if ~islogical(toPlot)
    % continuous variables
    switch numYvalues
        case 1                          % single vector (black)
            if isfield(method,'bar') && ~isempty(method.bar)
                % histogram
                bar(xValues, toPlot,'k')
                method.bar=[]; % avoid carry over between modules
            else
                % waveform
                plot(xValues, toPlot,'k')
            end
            xlim([0 method.objectDuration])
            if isfield(method,'zValuesRange') ...
                    && ~isempty(method.zValuesRange)
                ylim(method.zValuesRange)
                method.zValuesRange=[]; % avoid carry over between modules
            end
            if isfield(method,'yLabel') && ~isempty(method.yLabel)
                ylabel(method.yLabel, 'color', defaultTextColor)
                method.yLabel=[]; % avoid carry over between modules
            end
            
        case 2                          % 2 x N vector (black and red)
            plot(xValues, toPlot(1,:),'k'), % hold on
            plot(xValues, toPlot(2,:),'r'), % hold off
            xlim([0 method.objectDuration])
            if isfield(method,'zValuesRange') ...
                    && ~isempty(method.zValuesRange)
                ylim(method.zValuesRange)
                method.zValuesRange=[]; % avoid carry over between modules
            end
            if isfield(method,'yLabel')&& ~isempty(method.yLabel)
                ylabel(method.yLabel, 'color', defaultTextColor)
                method.yLabel=[]; % avoid carry over between modules
            end
            
        case 3                       % 3 x N vector (black red and green)
            % this is used for 1 channel DRNL output
            plot(xValues, toPlot(1,:),'k'), hold on
            plot(xValues, toPlot(2,:),'r'),  hold on
            plot(xValues, toPlot(3,:),'g'), hold off
            xlim([0 method.objectDuration])
            if isfield(method,'zValuesRange') ...
                    && ~isempty(method.zValuesRange)
                ylim(method.zValuesRange)
            end
            if isfield(method,'yLabel') &&  ~isempty(method.yLabel)
                ylabel(method.yLabel, 'color', defaultTextColor)
            end
            
        otherwise                       % >3 channels: surface plot
            % add  line to separate HSR and LSR
            if method.plotDivider && size(toPlot,1) > 2
                [r c]=size(toPlot);
                emptyLine=max(max(toPlot))*ones(2,c);
                halfway=round(r/2);
                toPlot=[toPlot(1:halfway,:); emptyLine; toPlot(halfway+1:end,:)];
            end
            
            % invert data for black on white matrix plotting
            if  method.blackOnWhite
                toPlot=-toPlot;
            end
            
            % matrix (analogue) plot
            if isfield(method,'forceLog') && ~isempty(method.forceLog)
                % positive values are put on log z-scale
                toPlot=toPlot+min(min(toPlot))+1;
                toPlot=log(toPlot);
                if isfield(method,'title')
                    method.title=[method.title '  (log scale)'];
                else
                    method.title= '(log scale)';
                end
            end
            
            %  zValuesRange
            if isfield(method,'zValuesRange') ...
                    && ~isempty(method.zValuesRange)
                clims=(method.zValuesRange);
                imagesc(xValues, yValues, toPlot, clims), axis xy; %NB assumes equally spaced y-values
            else
                % automatically scaled
                imagesc(xValues, yValues, toPlot), axis xy; %NB assumes equally spaced y-values
                
                if ~isfield(method,'zValuesRange')...
                        || isempty(method.zValuesRange)
                    method.zValuesRange=[-inf inf];
                end
                
                if method.blackOnWhite
                    % NB plotted values have negative sign for black on white
                    caxis([-method.zValuesRange(2) -method.zValuesRange(1)])
                else
                    caxis(method.zValuesRange)
                end
            end
            
            % xaxis
            % NB segmentation may shorten signal duration
            [r c]=size(toPlot);
            imageDuration=c*method.displaydt;
            xlim([0 imageDuration])
            
            % yaxis
            if isfield(method,'minyMaxy') && ~isempty(method.minyMaxy)
                ylim(method.minyMaxy)
            else
                if max(yValues)>min(yValues)
                    ylim([min(yValues) max(yValues)])
                end
            end
            
            % y-axis design yTickLabels
            if min(yValues)>1
                tickValues=[min(yValues) max(yValues)];
                tickLabels=num2str(tickValues');
                if method.plotDivider && size(toPlot,1) > 2
                    % show min/max yvalues with slight shift
                    yList=yValues;
                    yValues=1:length(yValues);
                    tickValues=[1 halfway-1 halfway+2 length(yValues)];
                    idx=[1 halfway halfway+1 length(yValues)];
                    tickLabels=num2str(yList(idx)');
                    imagesc(xValues, yValues, toPlot), axis xy;
                end
                
                set(gca,'ytick',tickValues)
                set(gca,'ytickLabel', strvcat(tickLabels))
                set(gca,'FontSize', method.defaultFontSize)
            end
            
    end
    
else	% is logical
    % logical implies spike array. Use raster plot
    [y,x]=find(toPlot);	%locate all spikes: y is fiber number ie row
    x=x*dt+method.timeStart;   % x is time
    plot(x,y, 'o', 'MarkerSize', rasterDotSize, 'color', 'k')
    if numYvalues>1
        set(gca,'yScale','linear')
        set(gca,'ytick', [1 numYvalues],'FontSize', method.defaultFontSize)
        % show lowest and highest BF value only
        set(gca,'ytickLabel', [min(yValues) max(yValues) ],'FontSize', method.defaultFontSize)
        if method.plotDivider
            % or use labels to identify fiber type
            set(gca,'ytickLabel', {'LSR', 'HSR'},'FontSize', method.defaultFontSize)
        end
        ylim([0 numYvalues+1])
    end
    xlim([0 method.objectDuration])
    if isfield(method,'yLabel') && ~isempty(method.yLabel)
        ylabel(method.yLabel,'FontSize', method.defaultFontSize, 'color', defaultTextColor)
    end
    
    % add  line to separate HSR and LSR
    if method.plotDivider
        [r c]=size(toPlot);
        halfWayUp=round(r/2);
        hold on
        plot([0 c*method.displaydt],[halfWayUp halfWayUp], 'b')
        hold off
    end
    
end

set(gca, 'xcolor', defaultAxesColor)
set(gca, 'ycolor', defaultAxesColor)

% add title
if isfield(method,'title') && ~isempty(method.title)
    title(method.title, 'FontSize', method.defaultFontSize, 'color', defaultTextColor)
end

% label axes
if ~isfield(method,'axes') || isempty(method.axes)
    % annotate the x-axis only if it is the last plot on a figure created by this utility
    set(gca,'xtick',[],'FontSize', method.defaultFontSize)
    if method.subPlotNo==method.numPlots
        if isfield(method,'xLabel') && ~isempty(method.xLabel)
            %               set(gca,'ActivePositionProperty','outerposition')
            %  xlabel(method.xLabel)
            xlabel(method.xLabel, 'FontSize', method.defaultFontSize, 'color', defaultTextColor)
        end
        set(gca,'xtickmode','auto') % add timescale to the lowest graph
    end
end

% add user labels to the y-axis if requested
if isfield(method,'yLabel') && ~isempty(method.yLabel)
    ylabel(method.yLabel, 'color', defaultTextColor)
end

% define color
if method.blackOnWhite, 	colormap bone, else 	colormap jet
end

% drawnow
