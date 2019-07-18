function exitPos = ultraBars(inData, colors, rowNames, ...
    horizonBars , nextAxis)

% Input:
% inData: a matrix and vector of plotting data
% colors: colors for each unique value of inData
% rowNames: a cell array for name for each row
% horizonBar: specifies whether to add a horizontal bar to each plot or not
 
% get the number of unique numbers and remove the 0 which means not plot
uniqueVars = unique(inData) ;
uniqueVars(uniqueVars == 0) = [] ;

% create the legend variable if the inData is a row vector
if size(inData,1) == 1
    lgdVar = split( num2str(uniqueVars) )';
end

% get the number of colours to plot
if ~exist('colors','var') || isempty(colors)
    colors = rand(length(uniqueVars), 3);
end

% check the are is a row name for each row in inData
if size(inData,1) > 1 % only valid for matrix data
    if size(inData,1) ~= size(rowNames,1)
        error('row names should contain a name for each row in plotting data')
    end
end

% check if the orizontal bars are required on the plot
if ~exist('horizonBars','var')
    horizonBars = false;
end

% check for color errors 
if length(uniqueVars) > size(colors,1) 
    error('A color must be specified for each value in Data')
end

% plot the figure for row vectors 
% I have separated the two peaces of code because one does not need a box
% around the data
if size(inData,1) == 1
    % figure()
    % axes('position',[0.1300,0.11,0.7,0.04]);
    for ii = 1:length(uniqueVars)
        % get the data for that values and plot it in that color
        plotData = double(ismember(inData,uniqueVars(ii) )) ;
        bar(plotData,'FaceColor',colors(ii,:),...
            'EdgeColor',[1 1 1] ,'BarWidth',0.9) ;
        hold on
    end
    set(gca,'GridColor',[1 1 1], 'XLim', [0.5 size(inData,2)+0.5 ], ...
        'XColor',[1 1 1] ,'YColor',[1 1 1],'XTick',[] , 'YTick',[],...
        'FontWeight','bold')
    % add a legend to the bar using the legendflex function
    
else
    % make a plot of multiple bar chart of top of each other
    % add the subtype plots to the heatmap using a loop
    % initialise some variables
    global plotTime
    figure(plotTime*100)
    % set(gcf,'position',[100,50,800,600])
        
    yInitial = 0.05; yPosSaved = yInitial;
    ySize = 0.02; increaseby = 0.022; % 0.44 0.4
    % change the size the bar if plotTime == 3
    if plotTime == 3
        ySize = 0.05; increaseby = 0.055;
    end
    xEndPos = 0.7 ; % 0.7750
    % loop over the rows and ascend by column
    for jj = 1:size(inData,1) % begin plotting from the bottem
        % define the next axis for the follow up plots
        if ~exist('nextAxis','var') || isempty(nextAxis)
            axes('position',[0.1300,yInitial,xEndPos,ySize]);
        elseif exist('nextAxis','var') && jj > 1
            axes('position',[0.1300,yInitial,xEndPos,ySize]);   
        else
            axes('position',nextAxis);
            yInitial = nextAxis(2) ;
            ySize = nextAxis(4) ; xEndPos = nextAxis(3) ;
            yPosSaved = yInitial ;
        end         
        for ii = 1:numel(uniqueVars) 
            plotData = double(ismember(inData(jj,:),uniqueVars(ii) )) ;
            bar(plotData,'FaceColor',colors(ii,:),'EdgeColor',[1 1 1] ,...
                'BarWidth',0.9) ;   
            hold on
            
            % add the name of the genes to the left of heatmap
            if exist('rowNames','var')
                dim = [0.02 yInitial 0.11 increaseby];
                annotation('textbox',dim,'String',rowNames{jj},...
                'FitBoxToText','on','FontSize',10,'EdgeColor','none',...
                    'HorizontalAlignment','right','FontWeight','bold',...
                    'VerticalAlignment','middle');
            end
        end
        % change the plot properties
        set(gca,'GridColor',[1 1 1], 'XLim', [0.5  size(inData,2)+0.5],...
            'XColor',[1 1 1] ,'YColor',[1 1 1],'YTickLabel',[],...
            'XTickLabel',[],'FontWeight','bold','YTick',[],'XTick',[])
        % increase the value to change the colors and plot positions
        yInitial = yInitial + increaseby;
    end
    % add a grey box to the plot usign the annotation function
    if plotTime ~= 3 % dont add the box if plotTime is == 3
    dim = [0.1300, yPosSaved, xEndPos, increaseby*size(inData,1)];
    annotation('rectangle',dim ,'Color',[0.5, 0.5, 0.5])
    end
    hold off 
   
    % plot the horizontal bars if they are required
    % prellocate the bar data size
    barhData = zeros(size(inData,1),numel(uniqueVars)) ;
    if horizonBars == true
        axes('position',[xEndPos+0.137, yPosSaved, 0.12, ...
            increaseby*size(inData,1)+0.001 ]);
        for kk = 1:numel(uniqueVars) 
            barhData(:,kk) = sum(inData == uniqueVars(kk),2) ;   
        end
        bar1 = barh(barhData,'stacked','BarWidth',0.85) ;
        % make sure there are no colors and spaces between the axis and the
        % first and last bar
        set(gca,'GridColor',[1 1 1], 'YLim', [0.5 size(inData,1)+0.52 ], ...
            'XColor',[0.5 0.5 0.5] ,'YColor',[1 1 1],'FontSize',10,...
            'YTick',[],'FontWeight','bold','Box','off','TickDir', 'out',...
            'LineWidth',1,'XAxisLocation','origin')
        
        % annoate the bar graph
        for ii = 1:size(barhData,2)
            set(bar1(ii),'FaceColor',colors(ii,:))
        end
        
    end
    exitPos = [0.1300,yInitial,xEndPos, ySize] ;
end

end