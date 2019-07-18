%% Pan-cancer Metabolic Gene Alterations - MAIN SCRIPT
% 
% MATLAB code to produce results present in the manuscript:Metabolic Gene 
% Alterations Across 32 Cancer Types and Their Clinical Significance

% change the directory and clear all the variable and the screen
% clc ; clear ; close all ;
% cd('/Users/sinkala/Documents/MATLAB/Metabolic Pan-Cancer')

% ================ Get Reactome Pathways Using an API ==================
% get the version of the database

fprintf('\n Getting Data from https://reactome.org/dev/content-service \n')

% Let’s start by retrieving the version of the database. It can be
% addressed by querying the “/data/database/version” method:
curl_get = sprintf('curl -X GET --header "Accept: text/plain" "https://reactome.org/ContentService/data/database/version"');
[~,databaseVersion] = system(sprintf('%s',curl_get));
fprintf('\n The Reactome Database Version is %d \n',...
    str2double(databaseVersion))

% Get all the stable IDs contained in the metabolic pathway reactome
% pathways
curl_get = sprintf('curl -X GET "https://reactome.org/ContentService/data/pathway/R-HSA-1430728/containedEvents/stId" -H "accept: text/plain"');
[~,reactomeIDs] = system(sprintf('%s',curl_get));
stableIDs = regexp(reactomeIDs, ',','split');
stableIDs = stableIDs' ;

%% Get All the Genes Involved in the Metabolic Pathways from Reactome  

fprintf('\n Getting Genes Involved in Metabolic Pathways \n')

% load the genes associated with each reactome pathway and change the name
% of the headers 
reactomePathways = readtable('ReactomePathways.txt');
reactomePathways.Properties.VariableNames = {'PathwayName','Genes'} ;

% load the text file that contains reactome pathways and only return the
% human pathway 
relations = readtable('ReactomePathwaysRelation.txt');
relations = relations( contains(relations.(1), 'R-HSA-'), :) ;

% get pathways that are associated with metabolism. Metabolic Pathways fall
% under the Top Level Pathway "Metabolism" = R-HSA-1430728 and print the
% number of metabolic processes associated with "Metabolism"
metabolicPathways = relations(contains(relations.(1),'R-HSA-1430728'),:) ;
metabolicPathways.Properties.VariableNames = {'Metabolism','subPathways'} ;
fprintf('The number of metabolic process is %d \n',height(metabolicPathways))

% get the pathways that are associated with each sub pathway
% preallocated the count for each process
appendCell = cell(height(metabolicPathways),4) ;
for ii = 1:height(metabolicPathways)
    % get the current sub-metabolic process: The first term is the
    % relations table where as the second term are the 16 metabolic
    % associated processes
    locProcess = contains(relations.(1), metabolicPathways.(2)(ii) ) ;
    appendCell(ii,2) = num2cell( sum(locProcess) ) ;
    
    % get the number of genes and genes names found in each of these
    % subprocesss and the name of the pathway
    curPathway = reactomePathways( contains( reactomePathways.(2), ...
        metabolicPathways.(2)(ii) ), : ) ;
    
    % add the name of the subpathway. followed the number of genes and
    % finally the actual gene names
    appendCell(ii,1) = strtrim( extractBefore( curPathway.(1), 'https') );
    
    % count the number of genes for that term and append the number to the
    % cell and also add the actuall gene names
    appendCell(ii,3) = num2cell( numel( split( extractAfter( ...
        curPathway.(2) ,metabolicPathways.(2)(ii) ) ) ) ) ;
    
    appendCell(ii,4) = strtrim( extractAfter(curPathway.(2) , ...
        metabolicPathways.(2)(ii) ) ) ;
end

metabolicPathways = addvars(metabolicPathways,appendCell(:,1),...
    appendCell(:,2),appendCell(:,3),appendCell(:,4),'NewVariableNames',...
    {'pathwayName','subSubPathway','countGenes','Genes'} );

%% Create a Table for all the Central Metabolic Pathways

fprintf('\n Creating Tables for Genes Involved in Metabolic Pathways \n')

% Get the subcategories of Carbohyrade Metabolism, Lipid Metabolism, Amino
% Acids Metabolism, TCA cycle and Integrated Metabolism 
centralCode = {'lipids',;'carbohydrates';'amino acid';'TCA';...
    'Integration'} ;
centralCode = metabolicPathways.subPathways(...
    contains(metabolicPathways.pathwayName, centralCode)  ) ;

centralPathways = relations(contains(relations.(1),centralCode ),:) ;
centralPathways.Properties.VariableNames = {'Pathway','subPathways'} ;
fprintf('The number of central metabolic process is %d \n',...
    height(centralPathways))

% get the pathways that are associated with each sub pathway
% preallocated the count for each process

appendCell = cell(height(centralPathways),4) ;
for ii = 1:height(centralPathways)
    % get the current sub-metabolic process: The first term is the
    % relations table where as the second term are the 16 metabolic
    % associated processes
    locProcess = contains(relations.(1), centralPathways.(2)(ii) ) ;
    appendCell(ii,2) = num2cell( sum(locProcess) ) ;
    
    % get the number of genes and genes names found in each of these
    % subprocesss and the name of the pathway
    curPathway = reactomePathways( contains(reactomePathways.(2), ...
        centralPathways.(2)(ii) ), : ) ;
    
    % add the name of the subpathway. followed the number of genes and
    % finally the actual gene names
    appendCell(ii,1) = strtrim( extractBefore( curPathway.(1), 'https') );
    
    % count the number of genes for that term and append the number to the
    % cell and also add the actuall gene names
    appendCell(ii,3) = num2cell( numel( split( extractAfter( ...
        curPathway.(2) ,centralPathways.(2)(ii) ) ) ) ) ;
    
    appendCell(ii,4) = strtrim( extractAfter(curPathway.(2) , ...
        centralPathways.(2)(ii) ) ) ;
end

centralPathways = addvars(centralPathways,appendCell(:,1),...
    appendCell(:,2),appendCell(:,3),appendCell(:,4), ...
    'NewVariableNames',{'pathwayName','subSubPathway','countGenes','Genes'});

% Finally add the Name of the pathway to the table 
appendCell = cell(height(centralPathways),1) ;
for ii = 1:length(centralCode)
    % get the current code and the current pathway name
    curName = metabolicPathways.pathwayName( ...
        contains(metabolicPathways.subPathways,centralCode(ii)) );
    
    % get the location of the central pathways with that current code
    locCentrals = contains(centralPathways.Pathway,centralCode(ii));
    
    %add the name to the cell array
    appendCell(locCentrals) = strtrim(curName);
end

centralPathways = addvars(centralPathways, appendCell,...
    'NewVariableNames','MainPathway','Before','Pathway') ;

%% ======== get all the genes involved in metabolic pathways ========

myGenes = [];
for ii = 1:height(metabolicPathways)
    curGenes = split(metabolicPathways.Genes(ii) );
    myGenes = [myGenes; curGenes] ;
end
% get the number of unique genes involved in metabolism and print to the
% screen 
myGenes = unique(myGenes) ;
fprintf('\nThe number of genes involved in metabolism = %d \n', ...
    length(myGenes) )

% clear some variables
clear appendCell centralCode curl_get curName curPathway ...
    databaseVersion ii locCentrals locProcess curGenes relations

%% ============== Get Data for My Genes from cBioPortal ==============

% Either load the cBioPortal data or get from cBioPortal If the data does
% not exist
if exist('metabolicData.mat','file')
    fprintf('\n Loading cBioPortal Data \n')
    load metabolicData.mat
    
    % change the CancerStudy to Upper
    mutations.CancerStudy = upper(mutations.CancerStudy);
    cnaData.CancerStudy = upper(cnaData.CancerStudy);
    clinicalData.CancerStudy = upper(extractBefore ( ...
        clinicalData.CancerStudy,'_TCGA')) ;
    % remove the duplicate sample Ids from the mutations data
    [~, theUnique] = unique(mutations.SampleIds);
    mutations = mutations(theUnique,:);
    
    % arrange the mutations data according to the cnaData
    [~, theLocs] = ismember(cnaData.SampleIds,mutations.SampleIds);
    mutations = mutations(theLocs, :);
    % throw in an assession
    assert(all(strcmp(mutations.SampleIds,cnaData.SampleIds)))
else
    fprintf(['\nGetting cBioPortal Data from http:/',...
        '/www.cbioportal.org/public-portal\n'])
    [mutations, cancerStudies, cnaData, mrna, clinicalData ] = ...
        getcBioPortalData(myGenes) ;
    
    % change the CancerStudy to Upper
    mutations.CancerStudy = upper(mutations.CancerStudy);
    cnaData.CancerStudy = upper(cnaData.CancerStudy);
%     clinicalData.CancerStudy = upper(extractBefore ( ...
%         clinicalData.CancerStudy,'_TCGA')) ;
    
    % save a copy of the data
    save('metabolicData.mat','mutations', 'cancerStudies', 'cnaData',...
        'mrna','clinicalData')
end

%% Perform Statistical Analyses and Plot Results

% Find the Number of participants for each cancer study
mutations.CancerStudy = categorical(mutations.CancerStudy);
participants = grpstats(mutations(:,1),'CancerStudy','count');
participants.Properties.VariableNames(1) = {'cancerTypeId'} ;
cancerStudies = innerjoin(cancerStudies,participants);
cancerStudies.mutationRate = cancerStudies.totalMutations./ ...
    cancerStudies.GroupCount ;

% Find the class of the cancer 
cancerClass = readtable('cancer classes.xlsx') ;
cancerClass.cancerTypeId = categorical(cancerClass.cancerTypeId);
cancerStudies = innerjoin(cancerStudies,cancerClass);
cancerStudies = movevars(cancerStudies,'cancerClass','before',...
    'totalMutations') ;
cancerStudies.cancerClass = categorical(cancerStudies.cancerClass);

clear cancerClass theUnique theLocs

%% get alterations of each metabolic pathways and put them in one table

% This is a table with columns as cancer types and row as metabolic
% pathways

% get the genes involved in a metabolic pathwy and return only these for
% the copy number data and mutations data
fprintf('\n Getting Pathway Alterations \n')
pathwayAlterations = findAlterationFreq(metabolicPathways,...
    mutations,cnaData) ;

% also get percentage of alteration involved in central metabolic pathways
centralPathwaysAlterations = findAlterationFreq(centralPathways,...
    mutations,cnaData) ;

%% plot the data for Figure 1

colorM  = cbrewer('seq', 'OrRd', 100); % BuPu

cgo = clustergram( pathwayAlterations{:,2:end}',...
    'columnlabels', cellstr(pathwayAlterations.CancerStudy),...
    'rowlabels', pathwayAlterations.Properties.VariableNames(2:end) ,...
    'colormap',redbluecmap ,'standardize','none',...
    'ColumnPDist','seuclidean');

% get a copy of the clustergram object and plot it for a sigle clustergram
cgoCopy = cgo ;
% get the colour bar on the plot 
cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
    'Color',{[1 1 0],[0.6 0.6 1]});
set(cgoCopy,'ColumnGroupMarker',cm ,'RowLabels',...
     strrep(pathwayAlterations.Properties.VariableNames(2:end)' ,'_','-'))

% make the figure very big to include the full row labels
cgAxes = plot(cgoCopy);
set(cgAxes, 'Clim', [0,100]) 

% get the low metabolism group: Get for new analysis!!!
lowMetabolism = clusterGroup(cgoCopy,30,'column','InfoOnly',true) ;
lowMetaBar = ~ismember(cgoCopy.ColumnLabels,lowMetabolism.ColumnNodeNames');

% add a heatmap on top to show the supervised classifier output
axes('position',[0.2423 0.730 0.5760 0.030] );  % 0.0250
ultraBars(double(lowMetaBar+1),[0 0.4509 0.7411; 0.8705882 0.490196 0],...
    {'Supertype'});
hold off


% ===== get the location of the data to produce the multiple plots ====
% there is a small bug where the row label seem to change by themselves 
set(cgo,'RowLabels',pathwayAlterations.Properties.VariableNames(2:end))
 
[~, locX] = ismember(flip(cgo.ColumnLabels), ...
    cellstr(pathwayAlterations.CancerStudy)) ;
[~, locY] = ismember(cgo.RowLabels, ...
    pathwayAlterations.Properties.VariableNames(2:end)) ;
heatData = pathwayAlterations{:,2:end}' ;
heatData = heatData(locY,locX);

% Produce Multiple Plots
figure(); clf
set(gcf,'position',[100,50,500,600]);
% the first number is how far the figure will be from the x-axis and the
% seceond number is now far the figure will be from the y-axis. The third
% number is was far the figure will run across the figure bar and the last
% number is far it will displaced allow the y-axis

axes('position',[0.20, 0.15, 0.74, 0.56]);
heatmap(cgo.ColumnLabels, flip(cgo.RowLabels), fliplr(flip(heatData)),...
    'Colormap',colorM,'ColorbarVisible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SET THE COLORS %%%%%%%%%%%%%%%%%%%%%%%%%%
rng(5) % 5 7
classesColors.cancerClass = ...
    rand(numel(unique(cancerStudies.cancerClass)),3) ;
classesColors.mutationRate = cool;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise the values for the for loopp
axes('position',[0.20, 0.715, 0.74,0.03]);
% change the arrangement in the cancer studies
[~, locX] = ismember(cgo.ColumnLabels, ...
    cellstr(cancerStudies.cancerTypeId)) ;
cancerStudies = cancerStudies(locX,:) ;
% now get the data of the cancer study
barData = cancerStudies.cancerClass' ;
heatmap(double(barData),'Colormap',classesColors.cancerClass,...
    'FontColor','none','ColorbarVisible','off');

axes('position',[0.2,0.75,0.74,0.12]);
barData = cancerStudies.mutationRate' ;
bar(barData ,'FaceColor',[0.5 .5 .5],'EdgeColor',[0.5 .5 .5]) ;
set(gca,'Box','off','XTick',[] ,'XLim',[0.5 32.5]);
ylabel('Mutation Rate')

clear locX locY barData colorM cgAxes cgoCopy

%% Correlation of Pathway Mutations to Patient Survival 

% get the low metabolism group: Get for new analysis!!!
lowMetabolism = clusterGroup(cgo,30,'column','InfoOnly',true) ;
lowMetabolismTumourNames = lowMetabolism.ColumnNodeNames';
lowMetaLoc = ~ismember(clinicalData.CancerStudy,...
    lowMetabolism.ColumnNodeNames') ;

% convert the logical into a cell and add it as a variable to the clinical
% data table 
metabolicStatus = ...
    strrep( (cellstr(num2str(lowMetaLoc))),'0','Low Metabolism');
metabolicStatus = strrep(metabolicStatus,'1','High Metabolism');
try % you dont want to get an error after the variable has been added to the table
    clinicalData = addvars(clinicalData,metabolicStatus,'Before',...
        'ABNORMAL_LYMPHOCYTE_PERCENT') ;
catch
end

% Get the Overall Survival Data and the Disease Free Survival Data and
% delete the missing low from the data
OsData = [clinicalData.OS_MONTHS, ...
    clinicalData.OS_STATUS ,num2cell(lowMetaLoc)] ;
OsData(any(cellfun(@isempty,OsData),2),:) = [] ;

fprintf('\nNumber of Low Metabolism Tumours = %d \n',...
    sum(cell2mat(OsData(:,3)) ))
fprintf('\nNumber of High Metabolism Tumours = %d \n',...
     sum(~cell2mat(OsData(:,3)) ))

% Anonotate the groups to either high metabolism or low metabolism and
% perferom K-M Survival Analysis 
groups = cellstr( num2str(cell2mat(OsData(:,3)) ) ) ;
groups = strrep(groups,'0','LM') ;
groups = strrep(groups,'1','HM') ;

% ============= PERFORM OVERALL SURVIVAL ANALYSIS =================
[~, ~, stats] = MatSurv( str2double(OsData(:,1)), ...
    lower(OsData(:,2)) , groups, 'NoRiskTable',true) ;

% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
    'TickDir','out')
ylevel = 0.65 ;
annotationData = str2double(OsData(:,1));
for jj = 1:numel(unique(groups))
    text(max(annotationData)/3,ylevel, ...
        sprintf('Median OS: (%s) = %g\n', ...
        stats.GroupNames{jj},stats.MedianSurvivalTime(jj)),...
        'FontWeight','bold','FontSize',10)
    ylevel = ylevel+0.04 ;
end  
title('\bf Overall Survival','FontSize',16)

% ============= PERFORM DISEASE FREE SURVIVAL ANALYSIS =================
% create for the data for disease free survival: first make the disease
% free data compatible with the matSurv function
dfsData = [clinicalData.DFS_MONTHS, ...
    strrep(regexprep(clinicalData.DFS_STATUS,'/(\w+)',''),...
    'Recurred', 'Relapsed'),num2cell(lowMetaLoc)] ;
dfsData(any(cellfun(@isempty,dfsData),2),:) = [] ;

% Anonotate the groups to either high metabolism or low metabolism and
% perferom K-M Survival Analysis 
groups = cellstr( num2str(cell2mat(dfsData(:,3)) ) ) ;
groups = strrep(groups,'0','LM') ;
groups = strrep(groups,'1','HM') ;

[~, ~, statsDFS] = MatSurv( str2double(dfsData(:,1)), ...
    dfsData(:,2), groups, 'NoRiskTable',true) ;

% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
    'TickDir','out')
ylevel = 0.65 ;
annotationData = str2double(dfsData(:,1));
for jj = 1:numel(unique(groups))
    text(max(annotationData)/3,ylevel, ...
        sprintf('Median DFS: (%s) = %g\n', ...
        statsDFS.GroupNames{jj},statsDFS.MedianSurvivalTime(jj)),...
        'FontWeight','bold','FontSize',10)
    ylevel = ylevel+0.04 ;
end  
title('\bf Disease Free Survival','FontSize',16)

% clear some variables 
clear ylevel annotationData groups OsData dfsData metabolicStatus ...
    lowMetaLoc lowMetabolism statsDFS stats ans jj heatData ...
    participants 

%% get Alteration of gene of Glycolysis and TCGA cycle.

% Munually curate these genes from the TCGA figures and the metabolic
% pahtways. These will be required to plot the pathways

% pathways are involve amino acids, lipids and carbohydrates
locCentral = ~cellfun(@isempty, regexp(centralPathways.MainPathway,...
    'lipids|amino acids|carbohydrates', 'match') ) ;
centrals = centralPathways.pathwayName(locCentral) ;
centrals2 = [true, ismember( ...
    centralPathwaysAlterations.Properties.VariableNames(2:end), ...
    matlab.lang.makeValidName(centrals)) ] ;
centralFinal = centralPathwaysAlterations(:,centrals2) ;

% arrange the centralFinal data according to the way the tumour types
% clustered in the first plot
[~,locCentral] = ismember(cgo.ColumnLabels ,centralFinal.CancerStudy);
centralFinal = centralFinal(locCentral,:) ;

% plot these similar to the first alterations figure
figure(); clf
set(gcf,'position',[100,50,700,500]);
axes('position',[0.20, 0.10, 0.74, 0.30]);
% get the length of lipids pathways. The variable names are taken for the
% original table and the curDataLoc will keep change per plot
% Produce Multiple Plots
curDataLoc = contains(centralPathways.MainPathway,'lipids') ;
curColor  = cbrewer('seq', 'Blues', 100); % Oranges Reds Greens
heatmap(centralFinal.CancerStudy, ...
    regexprep( centralPathways.pathwayName(curDataLoc), '\(+\w*)',''), ...
    centralFinal{:,2:sum(curDataLoc)+1}' ,...
    'Colormap',curColor,'ColorbarVisible','off');

% add the amino acids metabolism data to the plot and make sure you save
% the previous data length
axes('position',[0.20, 0.40, 0.74, 0.30]);
startPos = sum(curDataLoc)+ 2 ;
curDataLoc = contains(centralPathways.MainPathway,'amino acids') ;
endPos = startPos + sum(curDataLoc)-1 ;
curColor  = cbrewer('seq', 'Oranges', 100); % Oranges Reds Greens
h = heatmap(centralFinal.CancerStudy, ...
    regexprep( centralPathways.pathwayName(curDataLoc), '\(+\w*)',''), ...
    centralFinal{:,startPos:endPos}' ,...
    'Colormap',curColor,'ColorbarVisible','off');
hs = struct(h);
%now we can get at the private properties XAxis and YAxis
hs.XAxis.TickValues = [];

% add the carbohydrate metabolism data to the plot and make sure you save
% the previous data length
axes('position',[0.20, 0.70, 0.74, 0.25]);
startPos = endPos+1 ;
curDataLoc = contains(centralPathways.MainPathway,'carbohydrates') ;
endPos = startPos + sum(curDataLoc)-1 ;
curColor  = cbrewer('seq', 'Greens', 100); % Oranges Reds Greens
h = heatmap(centralFinal.CancerStudy, ...
    regexprep( centralPathways.pathwayName(curDataLoc), '\(+\w*)',''), ...
    centralFinal{:,startPos:endPos}' ,...
    'Colormap',curColor,'ColorbarVisible','off');
hs = struct(h);
%now we can get at the private properties XAxis and YAxis
hs.XAxis.TickValues = [];

%% Plot Glucose and B-oxidation Mutations 

fprintf('\n Creating Glycolysis, and Beta Oxidation Mutation Tables \n')

% get the reactome pathways involving Glycolysis, Gluconeogenesis, Beta
% Oxidation and Fatty Acid Biosynthesis
glucLipidCodes = {'R-HSA-70171','R-HSA-70263','R-HSA-75105','R-HSA-77289'};
glucLipidNames = {'Glycolysis',;'Gluconeogensis';...
    'LipidBiosythesis';'MitochondrialFattyAcidOxidation'} ;

% get the genes that are associated with each pathway
% preallocated the count for each process
for ii = 1:length(glucLipidCodes)
   curPathway = reactomePathways( ...
       contains(reactomePathways.Genes, glucLipidCodes{ii}), :) ;
   
   % return only the genes and get the alteration of these genes in each
   % cancer study and save the size of the geneAlterations table then also
   % delete the unmutated genes
   curGenes = split(strtrim( ...
       extractAfter(curPathway.Genes, glucLipidCodes{ii} )) ) ;
   geneAlters = mutations(:, ismember( ...
       mutations.Properties.VariableNames,['CancerStudy';curGenes] ) ) ;
   geneAlters(:, all(cellfun(@isempty ,geneAlters{:,2:end}) ,1) ) = [] ;
   
   % get the total alteration for each genes table and save them to a
   % strutured array
  glucLipidMutations.(glucLipidNames{ii}) = sum( ...
       ~cellfun(@isempty ,geneAlters{:,2:end}) , 1) / ...
       size(mutations,1) ;
   glucLipidMutations.(glucLipidNames{ii}) = ...
       [geneAlters.Properties.VariableNames(2:end)', ...
       num2cell(glucLipidMutations.(glucLipidNames{ii})*100 )' ] ;
   
   % make each study a row in the gene alteration table
   appendTable = cancerStudies(:,1) ;
   for jj = 1:height(cancerStudies)
       % get the alterations for the current cancer study and calculate the
       % pecentage of the alterations
       curStudy = geneAlters(ismember(geneAlters.CancerStudy,...
           cancerStudies.cancerTypeId(jj) ) ,:) ;
       
       % percentage of tumours that are altered for the cancer study
       perMut = array2table( ...
           sum(~cellfun(@isempty ,curStudy{:,2:end}) , 1) / ...
       cancerStudies.GroupCount(jj)*100 );
       appendTable(jj,2:width(curStudy) ) = perMut ; 
   end 
   
   % add variable names to the table
   appendTable.Properties.VariableNames(2:end) = ...
       curStudy.Properties.VariableNames(2:end) ;
  
   % delete the genes that are not mutated
   appendTable(:, [false , ~any(appendTable{:,2:end}) ] ) = [] ;
   
   % set the color
   if ii == 1 || ii == 2
       curColor  = cbrewer('seq', 'Greens', 100); % Oranges Reds Greens
   else
       curColor  = cbrewer('seq', 'Blues', 100); % Oranges Reds Greens
   end
   
   % plot the mutations on a chart
   % Produce Multiple Plots
   figure(); clf
%    set(gcf,'position',[100,100,1200,500]);
   % the first number is how far the figure will be from the x-axis and the
   % seceond number is now far the figure will be from the y-axis. The
   % third number is was far the figure will run across the figure bar and
   % the last number is far it will displaced allow the y-axis
   
   heatData = round( appendTable{:,2:end}' ,1); 
   heatData(heatData == 0 ) = NaN ;
   % axes('position',[0.1, 0.1, 0.8, 0.75]);
   if ii == 1
       h1 = heatmap(appendTable.cancerTypeId, ...
           appendTable.Properties.VariableNames(2:end),...
           heatData,'Colormap',curColor,'ColorbarVisible','off',...
           'MissingDataColor',[1 1 1]);
       h1.FontSize = 8;
       h1.Title = sprintf('%s Alterations',glucLipidNames{ii}) ;
   else
       h1 = heatmap(appendTable.cancerTypeId, ...
           appendTable.Properties.VariableNames(2:end),...
           heatData,'Colormap',curColor,'ColorbarVisible','off',...
           'MissingDataColor',[1 1 1]);
        h1.Title = sprintf('%s Alterations',glucLipidNames{ii}) ;
   end
  
end

clear heatData ii jj curColor perMut curStudy appendTable curPathway ...
    geneAlter curGenes centrals2 curDataLoc endPos glucLipidCodes ...
    glucLipidNames locCentral myGene MyGenes startPos ans h1

%% Plot the most mutated Gene Alterations 

% Who is the most mutated metabolic genes
fprintf('\n Finding the most mutated genes \n')
mostMutated = [mutations.Properties.VariableNames(3:end)' , ...
    num2cell( sum(~cellfun(@isempty,mutations{:,3:end}),1)' )] ;
mostMutated = sortrows(mostMutated,2 ,'descend') ;

% plot data of the most mutated gene
figure()
barh(categorical(mostMutated(1:30,1) ,flip( mostMutated(1:30,1) )), ...
    cell2mat( mostMutated(1:30,2) ), ...
    'FaceColor', [0.4666 0.1745 0.68823],'FaceAlpha',0.7) 
set(gca,'LineWidth',1,'FontSize',9,'Box','off','TickDir','in',...
    'FontWeight','bold')
xlabel('Genomic Alterations') ;
ylabel('Genes')
title('\bf Most Mutated Metabolic Genes','FontSize',12)

%% Overall Percentage Alterations for in All Tumours


% get the mutations and copy number data together
mutsAndCna = mutations;

% find the gene mutation frequencies
for ii = 3:width(mutsAndCna)
    % get the current genes 
    % check if that gene is present in the copy number data
    locGeneCNA = ismember(cnaData.Properties.VariableNames, ...
        mutsAndCna.Properties.VariableNames(ii) );
    
    % get the column of that genes
    if any(locGeneCNA)  
        curCol = cnaData{:,locGeneCNA} ~= 0 ;
        mutsAndCna.(ii)(curCol) = {'CNA'} ;
    end
end

% get the alterations for each genes 
eachCancerTotals = addvars(mutsAndCna(:,1) ,...
    sum(~cellfun(@isempty,mutsAndCna{:,3:end}),2),...
    'NewVariableNames','alterationRate') ;
eachCancerTotals.CancerStudy = categorical(eachCancerTotals.CancerStudy) ;
eachCancerTotals.alterProportion = logical(eachCancerTotals.alterationRate);

meanAltersAll = grpstats(eachCancerTotals,'CancerStudy') ;
meanAltersAll = sortrows(meanAltersAll,'mean_alterationRate','descend') ;
meanAltersAll.CancerStudy= reordercats( ...
    meanAltersAll.CancerStudy, cellstr(meanAltersAll.CancerStudy) ) ;

% for the average alterations of each pathways 
overallPercent = addvars(pathwayAlterations(:,1),...
    nanmedian(pathwayAlterations{:,2:end},2) ,'NewVariableNames',...
    'meanAlterationRate') ;
overallPercent = sortrows(overallPercent,'meanAlterationRate','descend');

%% produce Data to Use to Plot the Pathway

% also plot individuals genes on pathways using Draw.IO

% get the mutation frequencies to plot on the charts and remove the
% trailing spaces from the gene names
drawGenes = readtable('centralGenes.txt');
for ii = 1:width(drawGenes)
    drawGenes.(ii) = strip(drawGenes.(ii)) ;
end

if exist('drawAlterations.mat','file')
    fprintf('\n Loading mutations data for DrawIO \n')
    load('drawAlterations.mat')
else 
    fprintf(['\nGetting cBioPortal Data for DrawIO from http:/',...
        '/www.cbioportal.org/public-portal\n'])
    [drawMutations, ~, drawCNA, ~ ,~] = getcBioPortalData(drawGenes.Genes);
    
    % change the CancerStudy to Upper
    drawMutations.CancerStudy = upper(drawMutations.CancerStudy);
    drawCNA.CancerStudy = upper(drawCNA.CancerStudy);
    
    save('drawAlterations.mat', 'drawMutations', 'drawCNA')
end

% ===== try reading the table and see if i get cbioPortal results ====
drawCNA = readtable('centralCNA.txt');
drawMutations = readtable('centralMutations.txt');
drawMutations.Properties.VariableNames(2) = "SampleIds"  ;
drawCNA.Properties.VariableNames(2) = "SampleIds"  ;
% change the data type to double in the copy number data
for ii = 3:width(drawCNA)
    drawCNA.(ii) = str2double( drawCNA.(ii)) ;
    drawCNA.(ii)(drawCNA.(ii) >-2 & drawCNA.(ii) < 2 ) =  0 ;   
end
% replace the NA is the table
drawMutations = standardizeMissing(drawMutations,{'NA'});

% ==================================================================

% get the unique copy number mutations and copy number data
drawCNA = unique(drawCNA);
[~, locUnique] = unique(drawMutations.SampleIds);
drawMutations = drawMutations(locUnique,:) ;

% arrange the mutations and copy number data in the same way
[~, them] = ismember(drawCNA.Properties.VariableNames ,...
    drawMutations.Properties.VariableNames) ;
[~, them2] = ismember(drawCNA.SampleIds ,...
    drawMutations.SampleIds) ;
drawMutations = drawMutations(them2,them) ;

% find the gene mutation frequencies
for ii = 3:width(drawMutations)
    curCol = drawCNA.(ii) ~= 0 ;
    drawMutations.(ii)(curCol) = {'CNA'} ;
end
drawGeneAlters = sum(~cellfun(@isempty,drawMutations{:,3:end}),1) / ...
    height(drawMutations)*100;

% get the draws Genes table and add to it the frequencey of mutations for
% each gene
drawGeneAlters = table(drawCNA.Properties.VariableNames(3:end)',...
    drawGeneAlters','VariableNames',{'Genes','Proportions'} ) ;
drawGeneAlters = outerjoin(drawGenes(:,[1,3]),drawGeneAlters,...
    'MergeKey',true);

% finally get the percentage alteration of genes that belong to each class
drawGeneAlters.ProteinNames = categorical(drawGeneAlters.ProteinNames);
perDrawAlters = grpstats(drawGeneAlters(:,[1,3]),'ProteinNames','mean');
perDrawAlters.mean_Proportions = round(perDrawAlters.mean_Proportions,0);

% ============ Get the Alteration Frequencie in Each Cancer ==============
% This data will be use to plot a bar graph to put under the pathway

% clean up the names and remove the sample IDs
drawMutations.STUDY_ID = extractBefore(drawMutations.STUDY_ID ,'_');

% get the alterations for each genes 
drawBarData = addvars(drawMutations(:,1) ,...
    sum(~cellfun(@isempty,drawMutations{:,3:end}),2),...
    'NewVariableNames','alterationRate') ;
drawBarData.STUDY_ID = categorical(upper(drawBarData.STUDY_ID)) ;
drawBarData.alterDiscrete = logical(drawBarData.alterationRate);

meanAlters = grpstats(drawBarData,'STUDY_ID') ;
meanAlters = sortrows(meanAlters,'mean_alterDiscrete','descend') ;
meanAlters.STUDY_ID = reordercats( ...
    meanAlters.STUDY_ID, cellstr(meanAlters.STUDY_ID) ) ;

% plot the bar chart 
bar(meanAlters.STUDY_ID , meanAlters.mean_alterDiscrete*100,...
    'FaceColor',[0.5 .9 .5],'EdgeColor',[0.5 .9 .5],'FaceAlpha',0.8)
yt = get(gca, 'ytick');
ytl = strcat( strtrim( cellstr(num2str(yt'))),'%');
set(gca,'LineWidth',1,'FontSize',11,'Box','off','FontWeight','bold',...
    'XTickLabelRotation',90 ,'YTickLabel', ytl)
ylabel('Percentage Altered')
title('Alterations in Each Cancer Type','FontSize',15)

% add values to the chart 
y = meanAlters.mean_alterDiscrete*100 ;
x = 1:length(meanAlters.STUDY_ID);
for ii=1:numel(y)
    text(x(ii),y(ii), num2str(y(ii),'%0.0f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom')
end

% save the pentage of alterations to an excel file
% clear some of the variables
clear drawCNA drawMutations drawMutationsFreq drawCNAFreq them them2 ...
    drawGenes drawGeneAlters drawBarData drawMeanAlters meanAlter

%% Produce a Clustergram of Metabolic Alteration in All Tumour

% get all the clinical data related to the tumour stage
% tumourStages = [clinicalData(:,2), clinicalData(:,....
%     contains(clinicalData.Properties.VariableNames,'GRADE')) , ...
%     clinicalData(:, contains(clinicalData.Properties.VariableNames,...
%     'STAGE','IgnoreCase',true)) ] ;
% 
% rng(1)
% tumourStageColors = rand(20,3) ;
% 
% for ii = 1:height(cancerStudies)
%     curStudy = string(cancerStudies.cancerTypeId(ii)) ;
%     melanomaMuts = mutations(mutations.CancerStudy == curStudy,:) ;
%     melanomaCNA = cnaData(cnaData.CancerStudy == curStudy,:) ;
%     
%     melanomaGeneAlters =  findGeneAlterationFreq(metabolicPathways,...
%         melanomaMuts,melanomaCNA) ;
%     
%     cgoCancerStudies = clustergram( melanomaGeneAlters{:,2:end}',...
%         'columnlabels', cellstr(melanomaGeneAlters.SampleIds),...
%         'rowlabels', melanomaGeneAlters.Properties.VariableNames(2:end) ,...
%         'colormap',redbluecmap ,'standardize',1,...
%         'ColumnPDist','euclidean','Linkage','complete',...
%         'OptimalLeafOrder',1);
%     addTitle(cgoCancerStudies, sprintf('%s',curStudy),'FontSize',15)
%     
%     % ============== ADD the Tumour Grade to the Clustergram =============
%     % get the colour bar on the plot
%     cm = struct('GroupNumber',{1,2},'Annotation',{'Time1','Time2'},...
%         'Color',{[1 1 0],[0.6 0.6 1]});
% %     set(cgoCancerStudies,'ColumnGroupMarker',cm ,'RowLabels',...
% %         strrep(cgoCancerStudies.rowlabels,'_','-'))
%     
%     % make the figure very big to include the full row labels
%     cgAxes = plot(cgoCancerStudies);
%     
%     % get the tumour grade from the clinical data to use for the barplot
%     % for the current clustergram. It turns out that some the data dont not
%     % have clinical informations. Therefore I have to use an outerjoin
%     clustergramLabels = array2table(cgoCancerStudies.ColumnLabels', ...
%         'VariableNames',{'SampleIds'}) ;
%     curStages = tumourStages(ismember(tumourStages.SampleIds, ...
%         cgoCancerStudies.ColumnLabels') , :);
%     curStages = outerjoin(clustergramLabels,curStages ,...
%         'MergeKeys',true) ;
%     
%     % rearrange according to the way the samples clustered   
%     [~,reArrange] = ismember(cgoCancerStudies.ColumnLabels ,...
%         curStages.SampleIds);
%     curStages = curStages(reArrange,:);
%     
%     % throw in an assertions
%     assert( all(strcmp(curStages.SampleIds,cgoCancerStudies.ColumnLabels')));
%     
%     stageBar = curStages.GRADE ;
%     boxText = 'Grade' ;
%     
%     % if the cancer does not have the grade then get them from the other
%     % data
%     if all(cellfun(@isempty,stageBar))
%         stageBar = curStages.CLINICAL_STAGE ;
%         boxText = 'Clinical Stage' ;
%         if all(cellfun(@isempty,stageBar))
%             stageBar = curStages.AJCC_PATHOLOGIC_TUMOR_STAGE ;
%             boxText = 'AJCC Tumour Stage' ;
%         end
%         if all(cellfun(@isempty,stageBar))
%             stageBar = curStages.AJCC_CLINICAL_TUMOR_STAGE  ;
%         end
%         if all(cellfun(@isempty,stageBar))
%             stageBar = curStages.PATH_N_STAGE   ;
%             boxText = 'Pathology Stage' ;
%         end
%         if all(cellfun(@isempty,stageBar)) % if still nothing
%             fprintf('\n No Tumour Staging Data \n')
%             continue
%         end    
%     end
%     
%     % convert to categorical to allow ploting using alter bars and replace
%     % the missing categories (NaNs) with 0
%     stageBar = categorical(stageBar)' ;
%     % set the current colours 
%     curColors = tumourStageColors(1:length(categories(stageBar)),:);
%     stageBar = double(stageBar) ;
%     stageBar(isnan(stageBar)) = 0 ;
%     
%     % add a heatmap on top to show the supervised classifier output
%     axes('position',[0.2423 0.730 0.5760 0.030] );  % 0.0250
%     ultraBars(stageBar,curColors,{'Tumour Stage'} ,true);
%     
%     % ================ And the Names of the Staging ====================
%     % add the name of the cluster to the left of bar
%     dim = [0.1 0.730 0.11 0.030];
%     annotation('textbox',dim,'String', boxText,'FitBoxToText','on',...
%         'FontSize',12,'FontWeight','normal','EdgeColor','none',...
%         'HorizontalAlignment','right');
% 
%     % add the colours the left
%     hold on
%     if length(unique(stageBar)) < 6
%         xdistance = 0.03 ;
%         boxPos = [0.83 0.730 0.02 0.03] ;
%         boxColors = curColors;
%         stageBar(stageBar == 0) = [];
%         boxNumbers = unique(stageBar) ;
%         for jj = 1:length(unique(stageBar))
%             % add a box for each values 
%             annotation('rectangle',boxPos,'FaceColor',boxColors(jj,:), ...
%                 'EdgeColor',[1 1 1])
%             annotation('textbox',boxPos,'String',...
%                 num2str(boxNumbers(jj)),'FitBoxToText','on',...
%                 'FontSize',13,'FontWeight','bold','EdgeColor','none',...
%                 'HorizontalAlignment','center',...
%                 'VerticalAlignment','middle','Color',[1 1 1]) ;
%             boxPos(1) = boxPos(1) + xdistance ;
%         end
%     end
% end
% 
% clear stageBar curStages tumourStages boxColors boxNumbers ...
%     boxPos xdistance curColors boxText tumourStages tumourStageColors

%% Get Metabolic Alteration of ESCA and Check for Survival 

% ===== cluster the patients based on the mutational landscape =====
% get only the best cancer that show intratumour differences
% response data % T

oneCancer = 'ESCA' ;
melanomaMuts = mutations(mutations.CancerStudy == oneCancer,:) ; 
melanomaCNA = cnaData(cnaData.CancerStudy == oneCancer,:) ;

melanomaGeneAlters =  findGeneAlterationFreq(metabolicPathways,...
    melanomaMuts,melanomaCNA) ;

cgoMelanoma = clustergram( melanomaGeneAlters{:,2:end}',...
    'columnlabels', cellstr(melanomaGeneAlters.SampleIds),...
    'rowlabels', melanomaGeneAlters.Properties.VariableNames(2:end) ,...
    'colormap',redbluecmap ,'standardize',1,...
    'ColumnPDist','euclidean','Linkage','complete',...
    'OptimalLeafOrder',1);

cgAxes = plot(cgoMelanoma);
set(cgAxes, 'Clim', [-2, 2]) 
hold off

% ================ Create a Pathway for Abacavir Metabolism ========
% get the genes involved in abacavir metabolism

fprintf('\n Checking Mutations Different in %s Cancer \n',oneCancer)

% if the cancer is GBM get genes involved in inositol metabolism other wise
% get genes involved in Abacavir metabolism if the cancer is SKCM 

switch oneCancer
    case 'GBM'
        abcGenes = split( metabolicPathways.Genes(...
            contains( metabolicPathways.pathwayName,'Inositol')) ) ;
    case 'UCEC'
        abcGenes = split( metabolicPathways.Genes(...
            contains( metabolicPathways.pathwayName,'Inositol')) ) ;
    case 'CESC'
        abcGenes = split( metabolicPathways.Genes(...
            contains( metabolicPathways.pathwayName,'Inositol')) ) ;     
    otherwise
        abcGenes = split( metabolicPathways.Genes(...
            contains( metabolicPathways.pathwayName,'Abacavir')) ) ;
        
end

% get the high abacuvor metabolism group: Get for new analysis!!! .The
% abacavir groups has been saved as abacavir. Therefore, I load the data

% get the location of the names names that have mutations in more
% than 3% of all genes
InositolMutants = melanomaGeneAlters.SampleIds( ...
    melanomaGeneAlters.InositolPhosphateMetabolism > 5) ;
melanomaClinical = clinicalData(...
    clinicalData.CancerStudy == oneCancer,:);
locAbacavir =  ismember(melanomaClinical.SampleIds ,...
    InositolMutants) ;


% convert the logical into a cell and add it as a variable to the clinical
% data table 
metAbacavir = ...
    strrep( (cellstr(num2str(locAbacavir))),'1','Abacavir Pathway Mutated');
metAbacavir = strrep(metAbacavir,'0','No Mutations');

melanomaClinical = addvars(melanomaClinical,metAbacavir,'After',...
    'SampleIds','NewVariableNames','AbacavirMetabolism') ;

% ========================= Surival Analysis =====================
% survival analysis for abacavir metabolism vs non abacavir metabolism 

% Get the Overall Survival Data and the Disease Free Survival Data and
% delete the missing low from the data
OsData = [melanomaClinical.OS_MONTHS, ...
    melanomaClinical.OS_STATUS ,num2cell(locAbacavir)] ;
OsData(any(cellfun(@isempty,OsData),2),:) = [] ;

fprintf('\nNumber of Abacavir Unmutated Tumours = %d \n',...
    sum(cell2mat(OsData(:,3)) ))
fprintf('\nNumber of  Abacavir Mutant Tumours = %d \n',...
     sum(~cell2mat(OsData(:,3)) ))

% Anonotate the groups to either high metabolism or low metabolism and
% perferom K-M Survival Analysis 
groups = cellstr( num2str(cell2mat(OsData(:,3)) ) ) ;
groups = strrep(groups,'1','Abacavir Pathway Mutated') ;
groups = strrep(groups,'0','No Mutations') ;

% ============= PERFORM OVERALL SURVIVAL ANALYSIS =================
[~, ~, stats] = MatSurv( str2double(OsData(:,1)), ...
    lower(OsData(:,2)) , groups, 'NoRiskTable',true) ;

% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
    'TickDir','out')
ylevel = 0.65 ;
annotationData = str2double(OsData(:,1));
for jj = 1:numel(unique(groups))
    text(max(annotationData)/4,ylevel, ...
        sprintf('Median OS: (%s) = %g\n', ...
        stats.GroupNames{jj},stats.MedianSurvivalTime(jj)),...
        'FontWeight','bold','FontSize',12)
    ylevel = ylevel+0.04 ;
end  
title('\bf Overall Survival','FontSize',16)

% ============= PERFORM DISEASE FREE SURVIVAL ANALYSIS =================
% create for the data for disease free survival: first make the disease
% free data compatible with the matSurv function
dfsData = [melanomaClinical.DFS_MONTHS, ...
    strrep(regexprep(melanomaClinical.DFS_STATUS,'/(\w+)',''),...
    'Recurred', 'Relapsed'),num2cell(locAbacavir )] ;
dfsData(any(cellfun(@isempty,dfsData),2),:) = [] ;

% Anonotate the groups to either high metabolism or low metabolism and
% perferom K-M Survival Analysis 
groups = cellstr( num2str(cell2mat(dfsData(:,3)) ) ) ;
groups = strrep(groups,'1','Abacavir Pathway Mutated') ;
groups = strrep(groups,'0','No Mutations') ;

[~, ~, statsDFS] = MatSurv( str2double(dfsData(:,1)), ...
    dfsData(:,2), groups, 'NoRiskTable',true) ;

% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
    'TickDir','out')
ylevel = 0.65 ;
annotationData = str2double(dfsData(:,1));
for jj = 1:numel(unique(groups))
    text(max(annotationData)/4,ylevel, ...
        sprintf('Median DFS: (%s) = %g\n', ...
        statsDFS.GroupNames{jj},statsDFS.MedianSurvivalTime(jj)),...
        'FontWeight','bold','FontSize',12)
    ylevel = ylevel+0.04 ;
end  
title('\bf Disease Free Survival','FontSize',16)

clear metAbacavir locAbacavir aa abacavir abacavirNetwork ...
    ans annotationData abcMetabolism clustData colorM curStudy dfsData...
    fhsDFS fhOS geneAlters groups h heatData hs ii jj locX locY melanoma ...
    melanomaClinical melanomaCNA melanomaGeneAlters melanomaMuts ...
    myGenes OsData stats statsDFS ylevel

%% Load GDSC Data 

% ======================= GDSC Description =====================
% What genetic features are included in the ANOVA analysis for associations
% with drug response? 

% The ANOVA analysis currently correlates coding
% mutations and regions of recurrent copy number alteration with drug
% sensitivity data (IC50 values). These molecular alterations have been
% identified from the analysis of >11,000 patient tumours and subsequently
% mapped onto the cell lines for identifying molecular markers of drug
% sensitivity and resistance.

% ================== Load the Cell Line Mutation Data ===================
fprintf('\n Reading GDSC Data: This May Take While \n ')

fprintf('\n Reading GDSC Mutation Data\n ')
gdscMutations = readtable('GDSC Mutations.xlsx','Sheet',2);

fprintf('\n Reading GDSC Copy Number Data\n ')
gdscCNA = readtable('GDSC Copy Number Data.xlsx');
% clean up the gene symbols in the copy number data
gdscCNA.RegionIdentifier = strrep( erase( regexprep( extractAfter(...
    gdscCNA.RegionIdentifier, '(' ) ,'\,+\w*','') ,')') ,'-','_') ;
gdscCNA(cellfun(@isempty,gdscCNA.RegionIdentifier) ,:) = [] ;

fprintf('\n Reading GDSC Screened Compounds Data\n ')
gdscDrugs = readtable('GDSC_Screened_Compounds.xlsx');

fprintf('\n Reading GDSC Fitted Dose Response Curve Data\n ')
gdscDoseResponse = readtable('GDSC_fitted_dose_response.xlsx');

%% Check Drug Response of Abacavir Mutated Cell Line Vs Non-Mutants

fprintf('\n Checking InteraTumor Difference in Pathway \n')

% change the back to cancer types in the first column 
% get only the best cancer that show intratumour differences

% get the cancer study and return the data which has very high mutations
% rates
gdscStudies = unique(gdscMutations.CancerType) ;
gdscStudies(cellfun(@isempty,gdscStudies)) = [] ;
gdscStudiesBest = gdscStudies ;

%%%%%%%%%%%% THIS IS A HUGE LOOP FOR THE ALL THE CANCERS %%%%%%%%%%%%%%%%
for ii = 1 :height(metabolicPathways)
    fprintf('\n Finding Pathway Difference for %s',...
        metabolicPathways.pathwayName{ii}),
    
    for jj =  1 :length(gdscStudies)
        skcmMutation = gdscMutations(ismember(gdscMutations.CancerType, ...
            gdscStudies{jj}),:);
        skcmDoseRes = gdscDoseResponse( ...
            ismember(gdscDoseResponse.CELL_LINE_NAME, ...
            unique(skcmMutation.SAMPLE)),:);
        
        pathwayGenes = split( metabolicPathways.Genes(ii));
        
        skcmMutation.SAMPLE( ismember(skcmMutation.Gene, pathwayGenes) )  ;
        
        abcMutants = unique( ...
            skcmMutation.SAMPLE( ismember(skcmMutation.Gene,...
            pathwayGenes))) ;
        
        abcNotMutants = unique( ...
            skcmMutation.SAMPLE( ~ismember(skcmMutation.SAMPLE,...
            abcMutants))) ;
        
        % check if there is more than 4 mutatant. If not continue
        if numel(abcMutants) < 4 || numel(abcNotMutants) < 4
            fprintf('\n Less than 4 samples in on %s cancer subtypes\n',...
                gdscStudies{jj})
            continue
        end
        
        % add a column for the abcMutants in the table
        skcmDoseRes = addvars( skcmDoseRes ,...
            ismember(skcmDoseRes.CELL_LINE_NAME, abcMutants) ,...
            'After', 'COSMIC_ID','NewVariableNames','AbacavirMutant') ;
        
        % convert the drug names to categorical values
        skcmDoseRes.DRUG_NAME  = categorical(skcmDoseRes.DRUG_NAME);
        
        % perform t-test to compare the drug of responose of cells line
        % that have Abacavir pathway alterations vs those without pathway
        % alteration
        appendTable = unique(skcmDoseRes(:,[7,8])) ;
        toGoAppend = true(height(appendTable),1) ;
        for kk = 1:height(appendTable)
            % get the name of the drug
            curDrug = appendTable.DRUG_NAME(kk) ;
            
            % check that there are more than 2 cell lines being compared in
            % each analysis or comparison
            nMuts = numel(skcmDoseRes.LN_IC50(...
                skcmDoseRes.AbacavirMutant == true & ...
                skcmDoseRes.DRUG_NAME == curDrug)) ;
            nNonMuts = numel(skcmDoseRes.LN_IC50(...
                skcmDoseRes.AbacavirMutant  == false & ...
                skcmDoseRes.DRUG_NAME == curDrug));
            
            % continue the loop and set the rows to delete
            if nMuts < 3 || nNonMuts < 3
                fprintf('\n Less than 3 celllines in one group \n')
                toGoAppend(kk) = false ;
                continue
            end
            
            % perform a ttest with unequal variable assumed
            [p,h,stats] = ranksum( ...
                skcmDoseRes.LN_IC50(skcmDoseRes.AbacavirMutant == true ...
                & skcmDoseRes.DRUG_NAME == curDrug), ...
                skcmDoseRes.LN_IC50(skcmDoseRes.AbacavirMutant == false ...
                & skcmDoseRes.DRUG_NAME == curDrug) ,...
                'method','approximate') ;
            
            % get the means of the data
            meanValues = [ mean(  skcmDoseRes.LN_IC50(...
                skcmDoseRes.AbacavirMutant == true & ...
                skcmDoseRes.DRUG_NAME == curDrug) ),...
                mean( skcmDoseRes.LN_IC50(skcmDoseRes.AbacavirMutant ...
                == false & skcmDoseRes.DRUG_NAME == curDrug) ) ];
            
            % add to the table cant correct for multiple comparison if I
            % include the commented lines below
         
            appendTable(kk,3:8) = [ metabolicPathways.pathwayName(ii), ...
                gdscStudies(jj) ,...
                num2cell([meanValues, stats.zval ,p]) ] ;  
        end
        
        % delete the rows for which the ttest was not done
        appendTable = appendTable(toGoAppend,:) ;
        
        % add the variable names to the table
        appendTable.Properties.VariableNames(3:8) = ...
            {'MutatedPathway','CancerStudy','meanMuts','meanNonMuts',...
            'tValue','pValue'} ;
        
        % addd the p values to the table if only ther are some good results
        if ~all(isnan(appendTable.pValue))
            abacavirtTestResults = addvars( appendTable, ...
                mafdr(appendTable.pValue),'NewVariableNames','adjPvalue') ;
            abacavirtTestResults = sortrows(abacavirtTestResults,...
                'adjPvalue','ascend');
            
            % print the head of the data
            fprintf('\n Results for Study: %s\n',gdscStudies{jj} )
            head(abacavirtTestResults)
            
            % get the study with most significant changes in  the loop
            gdscStudiesBest(jj,ii+1) = ...
                num2cell(sum(abacavirtTestResults.adjPvalue < 0.05)) ;
            
            % save the significant results for the comparison between
            % tumours within a subtype with mutations and those without
            if ii == 1
            pathwayIntraTumourSig = ...
               abacavirtTestResults(abacavirtTestResults.pValue< 0.05,:);
            else 
                pathwayIntraTumourSig =...
                    vertcat(pathwayIntraTumourSig, ...
                    abacavirtTestResults( ...
                    abacavirtTestResults.pValue < 0.05,:) ) ;       
            end
            
        end
        
    end
    
end

writetable(pathwayIntraTumourSig, ...
    'pathway IntraTumour Drug Response.xlsx');
% convert to a table 
gdscStudiesBest = array2table(gdscStudiesBest ,'VariableNames',...
    ["cancerStudy", ...
    matlab.lang.makeValidName(metabolicPathways.pathwayName(...
    ~all(cellfun(@isempty,gdscStudiesBest))))' ] ) ;

clear appendTable ii p stats skcmDoseRes skcmMutation skcmCNA ...
    abcNotMutants abcMutants meanValues toGoAppend mMuts mNonMuts

%% Check Drug Response of High Metabolism Vs Low Metabolism

fprintf('\n Find Difference in Drug Action for Drugs \n')

% replace the / in the CancerType with nothing to make the names compareble
% to those by the TCGA
gdscMutations.CancerType = strrep(gdscMutations.CancerType,'/','');
DoseResponse = gdscDoseResponse ;

% get data that has only TCGA classified tumours from the GDSC
gdscOnlyTCGA = gdscMutations(ismember(gdscMutations.CancerType, ...
    cancerStudies.cancerTypeId),:);
onlyTCGASamples =unique( gdscMutations.SAMPLE(...
    ismember(gdscMutations.CancerType,cancerStudies.cancerTypeId),:) );
DoseResponse = DoseResponse( ...
    ismember(DoseResponse.CELL_LINE_NAME , onlyTCGASamples), :) ;

% return only genes involved in Low metabolism Genes
lowMetabolism = unique( ...
    gdscOnlyTCGA.SAMPLE( ismember(gdscOnlyTCGA.CancerType, ...
    lowMetabolismTumourNames) ) ) ;

% add a column for the abcMutants in the table 
DoseResponse = addvars( DoseResponse ,...
    ~ismember(DoseResponse.CELL_LINE_NAME, lowMetabolism) ,...
    'After', 'COSMIC_ID','NewVariableNames','Metabolism') ;

% convert the drug names to categorical values
DoseResponse.DRUG_NAME  = categorical(DoseResponse.DRUG_NAME);

% perform t-test to compare the drug of responose of cells line that have
% high metabolism and those with low metabolism
appendTable = unique(DoseResponse(:,[7,8])) ;
for ii = 1:height(appendTable)
    % get the name of the drug
    curDrug = appendTable.DRUG_NAME(ii) ;
    
    % perform a ttest with unequal variable assumed
    [~,p,~,stats] = ttest2( ...
        DoseResponse.LN_IC50(DoseResponse.Metabolism == true & ...
        DoseResponse.DRUG_NAME == curDrug), ...
        DoseResponse.LN_IC50(DoseResponse.Metabolism  == false & ...
        DoseResponse.DRUG_NAME == curDrug) ,...
        'Vartype','unequal') ;
    
    % get the mean values for the cell lines with metabolic alterations
    meanMetab = mean(...
        DoseResponse.LN_IC50(DoseResponse.Metabolism == true & ...
        DoseResponse.DRUG_NAME == curDrug)) ;
    
    meanNonMetab = mean(...
        DoseResponse.LN_IC50(DoseResponse.Metabolism  == false & ...
        DoseResponse.DRUG_NAME == curDrug)) ;
    
    % add to the table
    appendTable(ii,3:6) = num2cell([meanMetab,meanNonMetab,stats.tstat,p]);  
end

% add the variable names to the table 
appendTable.Properties.VariableNames(3:6) = ...
    {'meanHighMetabolism','meanLowMetabolism','tValue','pValue'} ;
MetabDrugResResults = addvars( appendTable, ...
    mafdr(appendTable.pValue,'BHFDR',true),...
    'After','pValue','NewVariableNames','FDR') ;
MetabDrugResResults = sortrows(MetabDrugResResults,'FDR','ascend');

% get the target pathways from the table of drugs
[~, themDrugs ] = unique(gdscDrugs.DRUG_NAME);
gdscDrugs = gdscDrugs(themDrugs,:);

% get the common drugs between the two table
gdscDrugs = gdscDrugs (...
    ismember(gdscDrugs.DRUG_NAME, MetabDrugResResults.DRUG_NAME) ,:) ;

% arrange both table in the same ways and then add a variable to the ttest
% results table
[~,arranGment] = ismember(...
    MetabDrugResResults.DRUG_NAME , gdscDrugs.DRUG_NAME);
gdscDrugs = gdscDrugs(arranGment,:) ;

% throw in a assertion 
assert(all(strcmp( cellstr(gdscDrugs.DRUG_NAME), ...
    cellstr(MetabDrugResResults.DRUG_NAME)) ))

% add the variable to the ttest results table and also get only the
% signficant results into one table
MetabDrugResResults = addvars(MetabDrugResResults , ...
    gdscDrugs.TARGET_PATHWAY, 'After','PUTATIVE_TARGET' ,...
    'NewVariableNames','TargetPathway') ;
sigDrugResults = MetabDrugResResults(MetabDrugResResults.FDR < 0.05, :);

% save the results to excel
writetable(MetabDrugResResults,'Metabolic Results.xlsx',...
    'Sheet','Drug Res- High vs Low Metab')
writetable(sigDrugResults,'Metabolic Results.xlsx',...
    'Sheet','Significant Results')

%% Plot three Box Plot for Metabolic Pathway A Targeting Drugs

% % get the the drug in the drug
metabolicDrugs = sigDrugResults(...
    contains( sigDrugResults.TargetPathway, 'Metabolism') & ...
    sigDrugResults.pValue < 0.05,: ) ;

% plot the graphs
for ii = 1:height(metabolicDrugs)
    % get the data and plot it on a graph
    figure
    plotData = DoseResponse.AUC(DoseResponse.DRUG_NAME == ...
        metabolicDrugs.DRUG_NAME(ii)) ;
    plotGroups = DoseResponse.Metabolism(DoseResponse.DRUG_NAME == ...
        metabolicDrugs.DRUG_NAME(ii) );
    
    % create the box plot
    boxplot(plotData, plotGroups) 
    % violinplot(plotData, plotGroups) 
    % scatterBoxPlots(plotData, plotGroups)
    
    % anonotate the figure
    set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
        'Box','off','XTickLabels',{'Low','High'})
    text(1.3, 0.95, sprintf('P = %0.4f', ...
        round(metabolicDrugs.pValue(ii),4) ) ,'FontSize',14)
    % add the ylabel 
    ylabel(sprintf('%s (%s) ',string(metabolicDrugs.DRUG_NAME(ii)) , ...
        metabolicDrugs.PUTATIVE_TARGET{ii} ))
    title('High Metabolism Vs Low Metabolism')
  
end

% IT TURN OUT THAT THE DATA DOES NOT LOOK ANY GOOD THEREFORE I FORGET ABOUT
% INLCUDING THESE DAT IN MY FINAL ANALYSIS 
clear metabolicDrugs plotData plotGroups ii gdscOnlyTCGA onlyTCGASamples

%% Get Drug response for each Drug that target Each pathways

fprintf('\n Find Difference in Drug Action for Pathways \n')

% add pathway information to the drug response data
for ii = 1:height(gdscDrugs)
    locDrug = ismember(DoseResponse.DRUG_NAME, gdscDrugs.DRUG_NAME(ii));
    DoseResponse.targetPathway(locDrug) = gdscDrugs.TARGET_PATHWAY(ii) ;
end

% convert to categorical so data acces is easy
DoseResponse.targetPathway = categorical(DoseResponse.targetPathway);

% perform t-test to compare the drug of responose of cells line that have
% high metabolic pathway alteration and those with low metabolic pathways
% alterations
appendTable = unique(DoseResponse(:,end)) ;
for ii = 1:height(appendTable)
    % get the name of the drug
    curPathway = appendTable.targetPathway(ii) ;
    
    % perform a ttest with unequal variable assumed
    [~,p,~,stats] = ttest2( ...
        DoseResponse.LN_IC50(DoseResponse.Metabolism == true & ...
        DoseResponse.targetPathway == curPathway), ...
        DoseResponse.LN_IC50(DoseResponse.Metabolism  == false & ...
        DoseResponse.targetPathway == curPathway) ,...
        'Vartype','unequal') ;
    
    % get the mean values for the cell lines with metabolic alterations
    meanMetab = mean(...
        DoseResponse.LN_IC50(DoseResponse.Metabolism == true & ...
        DoseResponse.targetPathway == curPathway)) ;
    
    meanNonMetab = mean(...
        DoseResponse.LN_IC50(DoseResponse.Metabolism  == false & ...
        DoseResponse.targetPathway == curPathway)) ;
    
    % add to the table
    appendTable(ii,2:5) = num2cell([meanMetab,meanNonMetab,stats.tstat,p]);  
end

% add the variable names to the table 
appendTable.Properties.VariableNames(2:5) = ...
    {'meanHighMetabolism','meanLowMetabolism','tValue','pValue'} ;

pathwaysDrugResResults = addvars( appendTable, ...
    mafdr(appendTable.pValue,'BHFDR',true),...
    'After','pValue','NewVariableNames','FDR') ;
pathwaysDrugResResults = sortrows(pathwaysDrugResResults,'FDR','ascend');

% save the data to excel
writetable(pathwaysDrugResResults ,'Metabolic Results.xlsx',...
    'Sheet','Pathway Drug Response Results')

% clear some of the variables
clear appendTable ii p stats skcmDoseRes skcmMutation skcmCNA ...
    abcNotMutants abcMutants themDrugs arranGment ...
    meanMetab meanNonMetab sigDrugResults curPathway qValue ...
    aa abcGene ans curDrug fhDFS locDrug lowMetabolism ...
    skcmDoseResponse smallResults stableIDs targetPathways them ...
    these

%% Produce some Code Plot for Signficant Pathways for Dose Response

% get the for most signifcation pathways between the high metabolism and
% low metabolism tumours 
mostSigPath = cellstr( pathwaysDrugResResults.targetPathway );
colNames = DoseResponse.Properties.VariableNames( ...
    contains(DoseResponse.Properties.VariableNames, ...
    {'Z_SCORE','targetPathway','CELL_LINE_NAME'} ) ) ;
selectData = DoseResponse( ismember( ...
    cellstr(DoseResponse.targetPathway) , mostSigPath), colNames) ;

% now get the data in a format that is usable by parrallelcoord plots
for ii = 1:length(mostSigPath)
    if ii == 1
        curPath = selectData(ismember(selectData.targetPathway, ...
            mostSigPath(ii)) , 1:2) ;
        curPath.Properties.VariableNames(2) = ...
            matlab.lang.makeValidName(mostSigPath(ii));
        % get the group stats
        curPath.CELL_LINE_NAME = categorical(curPath.CELL_LINE_NAME);
        curPath = grpstats(curPath(:,1:2),'CELL_LINE_NAME','nanmean') ;
        curPath = curPath(:,[1,3]) ;
        continue
    else
        curPath2 = selectData(ismember(selectData.targetPathway, ...
            mostSigPath(ii)) , 1:2) ;
        curPath2.Properties.VariableNames(2) = ...
            strrep( matlab.lang.makeValidName(mostSigPath(ii)) ,'_','');
        % get the group stats and then merge the tables
        curPath2.CELL_LINE_NAME = categorical(curPath2.CELL_LINE_NAME);
        curPath2 = grpstats(curPath2(:,1:2),'CELL_LINE_NAME','nanmean') ;
        % merge the table by using only the mean and actually cell line
        % names that are in column 1 and 3 of the two table
        curPath = outerjoin(curPath, curPath2(:,[1,3]),...
            'MergeKey',true);
    end
end

% add the metabolic status of the cell lines to be taken from the Drug
% respose data and join the two tables
metabClass = unique(DoseResponse(:,4:5));
metabClass.CELL_LINE_NAME = categorical(metabClass.CELL_LINE_NAME);
curPath = innerjoin(metabClass, curPath);

% produce a parrallel coodinate plots
figure
parallelcoords(curPath{:,3:end},'group',curPath.Metabolism,...
    'standardize','on','quantile',0.25,'LineWidth',1.5,...
    'labels', extractAfter( curPath.Properties.VariableNames(3:end) ,'_'))
set(gca,'YGrid','on','XGrid','on','XTickLabelRotation',45,...
    'LineWidth',1,'FontSize',10,'FontWeight','bold')

% ====== Produce the BAR graph that Darren suggests for Drug Action ======
forDrugsPlot = curPath(:,[1,3:end]) ;
drugBars = grpstats(forDrugsPlot ,'Metabolism',{'median','sem'}) ;

% get the data that contains nonmean 
meanVars = drugBars{ : , contains(drugBars.Properties.VariableNames,...
    'median')} ;
semVars = drugBars{ : , contains(drugBars.Properties.VariableNames,...
    'sem_nan')} ;

color = [0.00,0.45,0.74; 0.85,0.33,0.10];

figure
hold on
hb = bar(1:length(mostSigPath),meanVars') ;
% For each set of bars, find the centers of the bars, and write error bars
pause(0.1); %pause allows the figure to be created
for ib = 1:numel(hb)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    xData = hb(ib).XData+hb(ib).XOffset;
    er = errorbar(xData,meanVars(ib,:),semVars(ib,:),'k.');
    er.Color = color(ib,:) ; 
end
ylabel('Log IC50')
% add the tick and ticklabel
set(gca,'XTick',[1:length(mostSigPath)] ,'XTickLabel',mostSigPath ,...
    'XTickLabelRotation',45,'LineWidth',1)

% add the stars for pvalues
for ib = 1:numel(mostSigPath)
    %XData property is the tick labels/group centers; XOffset is the offset
    %of each distinct group
    % get the location of the first stars
    xVar = ib ;
    yVar = max(meanVars(:,ib)) + max(semVars(:,ib)) ;
    shift = 0.015 ;
    
    % change the shift if the values are both negative 
    if yVar < 0
       yVar = 0.02 ;
    end
    
    % get the p values and specificy the number of star to plot at that
    % points on the chart
    pvalue = pathwaysDrugResResults.pValue(ib) ;
    
    if pvalue < 0.001
        % plot three stars
        lineData = [xVar yVar; xVar yVar + shift; xVar yVar + 2*shift] ;     
    elseif pvalue < 0.01
         % plot two stars
         lineData = [xVar, yVar; xVar yVar + shift] ;
    elseif pvalue < 0.05
        lineData = [xVar, yVar];
        % plot two stars
    else
        % plot nothing
        continue;
    end
    
    % finally produce the line plot
    line(lineData(:,1),lineData(:,2),'LineStyle','none',...
        'Marker','pentagram','MarkerSize',5,'Color','k',...
        'MarkerFaceColor','k')
end
hold off

% Clear some Variables
clear forDrugsPlot hb ib xData semVar meanVars drugBars curPath ...
    colNames mostSigPath curPath2 selectData metabClass ...
    themMetab curCol ii locMetab locDrug color pvalue xyVar ...
    shift xVar yVar ib
%% Get Metabolic Pathway Alterations in GDSC Cell Lines

% By checking the mutations landscape of the cell line verses the drug
% response profile of the cell lines

try
    % save a copy of the TCGA cell line classification
    tcgaDescript = unique(gdscMutations(:,[1,3]));
    tcgaDescript.Properties.VariableNames(1) = "SampleIDs" ;
    
    % get the mutations and copy number data into tables
    fprintf('\n Processing GDSC Cell Lines mutation Data \n')
    gdscMutations = processMAF(gdscMutations) ;
    fprintf('\n Processing GDSC Cell Lines Copy Number Data \n')
    gdscCNA = processMAF(gdscCNA);
    
    % join the table with the tcga descriptions 
    gdscMutations = outerjoin(tcgaDescript, gdscMutations, 'MergeKey',true);
    gdscCNA = outerjoin(tcgaDescript, gdscCNA, 'MergeKey',true);
    
    % delete the row which do not have a TCGA cancer type
    gdscMutations(cellfun(@isempty ,gdscMutations.CancerType),:) = [];
    gdscCNA(cellfun(@isempty ,gdscCNA.CancerType),:) = [];
    
    % empty columns are causing the error: fix tomorrow
    gdscMutations(:, all(cellfun(@isempty,gdscMutations{:,3:end}))) = [];
    gdscCNA(:, all(cellfun(@isempty ,gdscCNA{:,3:end})) ) = [];
    
    % move the cancerType variable to the first column and rename the
    % variable to cancerstudy
    gdscMutations.Properties.VariableNames(2) = "CancerStudy" ;
    gdscMutations = movevars(gdscMutations,'CancerStudy',...
        'Before','SampleIDs');
    gdscCNA.Properties.VariableNames(2) = "CancerStudy" ;
    gdscCNA=movevars(gdscCNA,'CancerStudy','Before','SampleIDs');  
    
    % convert the copy number data to either -2 and 2 for amplification and
    % deletions
    for jj = 3:width(gdscCNA)
        gdscCNA.(jj)(strcmpi(gdscCNA.(jj),'Deletion')) = ...
            num2cell(-2) ;
        gdscCNA.(jj)(strcmpi(gdscCNA.(jj),'Amplification')) = ...
            num2cell(2) ;
    end

catch
end

% get the pathway alterations and arrange them according the clustering of
% the data
fprintf('\n Getting Pathway Alterations in GDSC Cell Lines \n')
pathwayAltersGDSC = findAlterationFreq(metabolicPathways,...
    gdscMutations,gdscCNA ) ;

% ================= plot the data for Figure 1 ======================
colorM  = cbrewer('seq', 'OrRd', 100); % BuPu
cgo2 = clustergram( pathwayAltersGDSC{:,2:end}',...
    'columnlabels', cellstr(pathwayAltersGDSC.CancerStudy),...
    'rowlabels', pathwayAltersGDSC.Properties.VariableNames(2:end) ,...
    'colormap',redbluecmap ,'standardize','none',...
    'ColumnPDist','seuclidean');

cgAxes = plot(cgo2);
set(cgAxes, 'Clim', [0,100]) 
hold off

[~, locX] = ismember(flip(cgo2.ColumnLabels), ...
    cellstr(pathwayAltersGDSC.CancerStudy)) ;
[~, locY] = ismember(cgo2.RowLabels, ...
    pathwayAltersGDSC.Properties.VariableNames(2:end)) ;
heatData = pathwayAltersGDSC{:,2:end}' ;
heatData = heatData(locY,locX);

figure()
% the first number is how far the figure will be from the x-axis and the
% seceond number is now far the figure will be from the y-axis. The third
% number is was far the figure will run across the figure bar and the last
% number is far it will displaced allow the y-axis

heatmap(cgo2.ColumnLabels, flip(cgo2.RowLabels), fliplr(flip(heatData)),...
    'Colormap',colorM,'ColorbarVisible','off');

%% Plot the Cancer Cell lines and Primary Tumours Mutations Side by Side

% get the intersection between GDSC cell line class and TCGA 
commonCancer = intersect(pathwayAlterations.CancerStudy ,...
    pathwayAltersGDSC.CancerStudy);
commonGDSC = pathwayAltersGDSC( ...
    ismember(pathwayAltersGDSC.CancerStudy, commonCancer),:) ;
commonTCGA = pathwayAlterations( ...
    ismember(pathwayAlterations.CancerStudy, commonCancer),:) ;

% change the cancer names in GDSC pathway alterations to add cell lines
commonGDSC.CancerStudy = strcat(cellstr(commonGDSC.CancerStudy),'cl') ;

% combine the two tables 
TCGA_GDSCpathways = vertcat(commonTCGA, commonGDSC) ;
TCGA_GDSCpathways.CancerStudy = cellstr(TCGA_GDSCpathways.CancerStudy);
TCGA_GDSCpathways = sortrows(TCGA_GDSCpathways,'CancerStudy','ascend') ;

% re-arrange the column according the arrange in the clustered data
[~, locY] = ismember(['CancerStudy';cgo2.RowLabels], ...
   TCGA_GDSCpathways.Properties.VariableNames) ;
TCGA_GDSCpathways = TCGA_GDSCpathways(:,locY);

% finally create a heat map 
figure; clf
axes('position',[0.20, 0.10, 0.65, 0.75]);
heatmap(TCGA_GDSCpathways.CancerStudy,...
    flip(TCGA_GDSCpathways.Properties.VariableNames(2:end)),...
    flip(TCGA_GDSCpathways{:,2:end}'),...
    'Colormap',colorM,'ColorbarVisible','off');
% add a bargrap to the plot of the mutation rates in each tumours calss
axes('position',[0.20, 0.855, 0.65, 0.1]);
% colorTCGA = [0  0.4470 0.7410; 0.8500 0.3250 0.0980 ];
colorTCGA = [0 0.7509 0.9411; 0.905882 0.790196 0.2];
barData = mean(TCGA_GDSCpathways{:,2:end},2)' ;
b1 = bar(barData ,'FaceColor','flat','FaceAlpha',0.7) ;
set(gca,'Box','off','XTick',[] ,'XLim',[0.5 length(barData)+0.5] ,...
    'YTick',[30 60],'YTickLabel',{'30%','60%'});
for ii = 1:length(barData)
    if rem(ii,2) ~= 0 % if the number is not even
        b1.CData(ii,:) = colorTCGA(1,:);
    else
        b1.CData(ii,:) = colorTCGA(2,:);
    end
end

% add a bar graph to the left of hte plot which show the percentage of
% mutations in each pathway for each tumour class
axes('position',[0.853, 0.10, 0.06, 0.75]);
gdscRows = contains(TCGA_GDSCpathways.CancerStudy,'cl') ;
barData = [ mean(TCGA_GDSCpathways{ ~gdscRows,2:end},1)'/2, ...
     mean(TCGA_GDSCpathways{gdscRows,2:end},1)'/2 ] ;
b2 = barh(barData,'stacked') ;
b2(1).FaceColor = colorTCGA(1,:) ; b2(1).FaceAlpha = 0.7;
b2(2).FaceColor = colorTCGA(2,:) ; b2(2).FaceAlpha = 0.7;
set(gca,'Box','off','YTick',[] ,'YLim',[0.5 length(barData)+0.5] ,...
    'XTick',[50 100],'XTickLabel',{'50%','100%'});

% clear some variables
clear tcgaDescript locX locY barData cgAxes colorTCGA commonGDSC...
     b1 barData colorTCGA locUnique ...
    gdscRows mostSigPath patients

% ======== Produce Some Bar Graphs for Significant Pathway ============

% some drug have been shown to be more active on cell lines with lower
% metabolic pathway alteration than those of tumours with high metabolic
% pathways alterations. Here I produce plots for all the drugs that belong
% to a specific family comparing the drug action across cancers that belong
% to each of the two metabolic pathways phenotypes.

% add CGA classification to the drug reponse data
[these, them] = ismember(gdscDoseResponse.CELL_LINE_NAME,...
    gdscMutations.SampleIDs);
% delete the zeroes from the index values 
them = them(these) ;
gdscPathwayResponse = gdscDoseResponse(these,:) ;
gdscPathwayResponse = addvars(gdscPathwayResponse ,...
    gdscMutations.CancerStudy(them),'After', 'CELL_LINE_NAME',...
    'NewVariableNames','CancerType');

% add the target pathways to the gdscPathwayResponse data
[these, them] = ismember(gdscPathwayResponse.DRUG_NAME,...
    gdscDrugs.DRUG_NAME);
% delete the zeroes from the index values 
them = them(these) ;
gdscPathwayResponse = gdscPathwayResponse(these,:) ;
gdscPathwayResponse = addvars(gdscPathwayResponse ,...
    gdscDrugs.TARGET_PATHWAY(them),'After', 'PUTATIVE_TARGET',...
    'NewVariableNames','TargetPathway');

% add the low metabolism tag to the cell lines in the table 
these = ismember(gdscPathwayResponse.CancerType,lowMetabolismTumourNames);
% delete the zeroes from the index values 
gdscPathwayResponse = addvars(gdscPathwayResponse ,...
    ~these,'After','CancerType','NewVariableNames','MetabolicStatus');

% Finally return only the cell lines that are members of the TCGA
% classification scheme
gdscPathwayResponse = gdscPathwayResponse( ...
    ismember(gdscPathwayResponse.CancerType,cellstr(commonCancer)) ,: );

%% ================ Produce the Plots in a Loop =========================

targetPathways = unique(gdscPathwayResponse.TargetPathway) ;
for ii = 1:length(targetPathways)
    % get the Cancer Types Metabolic Status and IC50 values for the current
    % pathway
    curTable = gdscPathwayResponse ( ...
        ismember(gdscPathwayResponse.TargetPathway,...
        targetPathways(ii) ) , [5,6,13] );
    curTable = sortrows(curTable,'MetabolicStatus','descend') ;
    
    % set the color of the box plots
    % first get the number of low metabolism tumours
    colorCutoff = ...
        sortrows( unique(curTable(:,[1,2])),'MetabolicStatus','ascend') ;
    color = [0 0.4509 0.7411; 0.8705882 0.490196 0];
    
    figure(100+ii)
    axes('position',[0.05, 0.1, 0.70, 0.85]);
    boxplot(curTable.LN_IC50, curTable.CancerType,'ColorGroup',...
        curTable.MetabolicStatus);
    % change the figure style
    figure(100+ii)
    ylabel('Ln IC50','FontWeight','bold')
    xlabel('Cancer Type','FontWeight','bold')
    title(sprintf('Dose Response of %s Targeting Drugs by Cancer Type',...
        targetPathways{ii}),'FontSize',14','FontWeight','bold')
    set(gca,'Box','off','LineWidth',1.5);
    h = findobj(gca,'Tag','Box') ;
    for jj = 1:length(h)
        switch colorCutoff.MetabolicStatus(jj)
            case false
                patch(get(h(jj),'XData'),get(h(jj),'YData'),...
                    color(1,:),'FaceAlpha',.5,'LineWidth',1);
            case true
                patch(get(h(jj),'XData'),get(h(jj),'YData'),...
                    color(2,:),'FaceAlpha',.6,'LineWidth',1);
        end
    end
    
    % plot a box plots
    axes('position',[0.8, 0.1, 0.18, 0.85]);
    boxplot(curTable.LN_IC50, curTable.MetabolicStatus)
    % color the box plots 
    color = [0.8705882 0.490196 0 ; 0 0.4509 0.7411];
    h4 = findobj(gca,'Tag','Box') ;
    for j=1:length(h4)
        patch(get(h4(j),'XData'),get(h4(j),'YData'),color(j,:),...
            'FaceAlpha',.6,'LineWidth',1);
    end
    % change the  box properties and add the pvalue to the plot for the
    % comprision 
    set(gca,'Box','off','LineWidth',1.5, 'XTick',[1 2],...
        'XTickLabel',{'Low Mutations','High Mutations'});
    title('Overall Group Comparision','FontSize',14','FontWeight','bold')
    % get the p value
    curPvalue = pathwaysDrugResResults.pValue( ...
        ismember(pathwaysDrugResResults.targetPathway, ...
        targetPathways(ii)) ) ;
    text(1.35,max(curTable.LN_IC50),sprintf('p = %0.4f',curPvalue))
end

clear h ii jj color curTable colorCutoff

%% ============= Produce the Plots For Signficant Drugs ===================

% prelloacte the drug table
appendTable = gdscDrugs(:,[2,4:end]) ;
drugs = gdscDrugs.DRUG_NAME ;
for ii = 1:length(drugs)
    % get the Cancer Types Metabolic Status and IC50 values for the current
    % pathway
    curTable = gdscPathwayResponse( ...
        ismember(gdscPathwayResponse.DRUG_NAME,...
        drugs(ii) ), [5,6,13] ) ;
    curTable = sortrows(curTable,'MetabolicStatus','descend') ;
    curResponse = curTable.LN_IC50 ;
    
    % perform a t-test if the results are statisitically significant then
    % plot the data alse continue
    [~,curPvalue,~,stats] = ...
        ttest2(curResponse(curTable.MetabolicStatus == 1),...
        curResponse(curTable.MetabolicStatus == 0),...
        'Vartype','unequal') ;
    
    fprintf('\nThe Current P Values for Drug %s is %d\n',drugs{ii},curPvalue)
    % create a table of the drug response to the cell lines to the
    % anticancer drugs
    meanMetabHigh = mean(curResponse(curTable.MetabolicStatus == 1)) ; 
    meanMetabLow = mean(curResponse(curTable.MetabolicStatus == 0)) ;
    nMetabHigh = length(curResponse(curTable.MetabolicStatus == 1));
    nMetabLow = length(curResponse(curTable.MetabolicStatus == 0));
    
    appendTable(ii,[4:9]) = num2cell([meanMetabHigh, meanMetabLow,...
       nMetabHigh,  nMetabLow,stats.tstat,curPvalue]);
    
    if curPvalue > 0.0001
        clear curPvalue
        fprintf('\n Loop exist because of high p value\n');
        continue
    end
    
    figure()
    axes('position',[0.05, 0.1, 0.70, 0.85]);
    hold on
    boxplot(curResponse, curTable.CancerType,'ColorGroup',...
        curTable.MetabolicStatus)
    % set the color of the box plots
    % first get the number of low metabolism tumours 
    rng(1,'twister') % for color reproducibility
    color = rand(height(curTable),3) ;
    h = findobj(gca,'Tag','Box');
    for jj = 1:length(h)
        patch(get(h(jj),'XData'),get(h(jj),'YData'),...
            color(jj,:),'FaceAlpha',.8,'LineWidth',1);
    end
    set(gca,'Box','off','LineWidth',1.5);
    ylabel('Ln IC50','FontWeight','bold')
    xlabel('Cancer Type','FontWeight','bold')
    title(sprintf('IC50 Distriubtion for %s by Cancer Type',...
        drugs{ii}),'FontSize',14','FontWeight','bold')
    axes('position',[0.8, 0.1, 0.18, 0.85]);
    boxplot(curResponse, curTable.MetabolicStatus)
    % color the box plots 
    color = [0.8705882 0.490196 0 ; 0 0.4509 0.7411];
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),...
            'FaceAlpha',.7,'LineWidth',1);
    end
    % change the  box properties and add the pvalue to the plot for the
    % comprision 
    set(gca,'Box','off','LineWidth',1.5, 'XTick',[1 2],...
        'XTickLabel',{'Low Mutations','High Mutations'});
    title('Overall Group Comparision','FontSize',14','FontWeight','bold')
    % get the p value
    text(1.35,max(curResponse),sprintf('p = %d',curPvalue))
    hold off
end
DrugDrugResResults = appendTable ;

% add the variable names
DrugDrugResResults.Properties.VariableNames(4:9) = ...
    {'meanHighMetabolism','meanLowMetabolism','nHigh','nLow',...
    'tValue','pValue'} ;

% add the adjusted p values to the table 
DrugDrugResResults = addvars(DrugDrugResResults, ...
    mafdr(DrugDrugResResults.pValue,'BHFDR',true),...
    'After','pValue','NewVariableNames','FDR') ;
DrugDrugResResults  = sortrows(DrugDrugResResults ,'FDR','ascend');

sigDrugResults = DrugDrugResResults(DrugDrugResResults.FDR < 0.05, :);

% save the results to excel
writetable(DrugDrugResResults,'Metabolic Results.xlsx',...
    'Sheet','Drug Res- High vs Low Metab')
writetable(sigDrugResults,'Metabolic Results.xlsx',...
    'Sheet','Significant Results')

% this assertion comfirms that the two dataset are exactly the same 
assert(all(gdscPathwayResponse.MetabolicStatus == DoseResponse.Metabolism))

% clear some of the variables
clear h ii jj color curTable colorCutoff appendTable meanMetabHigh...
    meanMetabLow curPvalue stats

%% mRNA based clustering the Two Cancer Groups

fprintf('\n Producing tnse Plots for the mRNA data \n')

% plot some tnse with the data later plotted in tableau of the tsne similar
% to what I produced for the clinical trials
mrna.CancerStudy = upper(mrna.CancerStudy) ;

try % try incase this has already happened 
    mrna = addvars(mrna,...
        ~ismember(mrna.CancerStudy,lowMetabolismTumourNames),...
        'After','CancerStudy','NewVariableNames','MetabolicStatus');
catch
end

% convert all the variable in the table to double
toGo = false(1,width(mrna));
for ii = 4:width(mrna)
    if iscell(mrna.(ii)) % if its a cell array
        if isempty( cell2mat(mrna.(ii) ) ) %if here is nothing in the cell
            toGo(ii) = true ;
        else 
            mrna.(ii) = cell2mat(mrna.(ii)) ;
        end
    end
end
% now delete the to go rows
mrna(:,toGo) = [] ;

mrna4Corr = mrna(:,[1,2,4:end]);

% obtain expression measurements
mrnaExpr = mrna{:,4:end}';
genes = mrna.Properties.VariableNames(4:end) ;
groups = mrna.MetabolicStatus ;
mrna.CancerStudy = categorical(mrna.CancerStudy);
groupsCancer = mrna.CancerStudy;

% remove nan values
nanIndices = any(isnan(mrnaExpr),2);
mrnaExpr(nanIndices,:) = [];
genes(nanIndices) = [];
numel(genes)

% ======================== Filter Out Genes =====================

% Gene profiling experiments typically include genes that exhibit little
% variation in their profile and are generally not of interest. These genes
% are commonly removed from the data.
for times = 1 % 5 4,3,2 also does good plot but one groups is too separated
    mask = genevarfilter(mrnaExpr);
    
    mrnaExpr = mrnaExpr(mask,:);
    genes = genes(mask);
    numel(genes)
    
    % filter out genes below a certain fold change threshold
    [~,mrnaExpr,genes] = ...
        genelowvalfilter(mrnaExpr,genes,'absval',log2(3));
    numel(genes)
    
    % filter genes below a certain percentile: VERY POWERFUL discriminant
    [~,mrnaExpr,genes] = ...
        geneentropyfilter(mrnaExpr,genes,'prctile',40);
    numel(genes)
end

% ==================================================================== % 
fprintf('\n Now Running tSNE for 2 Dimesions \n')
% use maybe spectrial clustering here
rng default % for reproducibility
yValues = tsne(mrnaExpr','Algorithm','barneshut','Distance','spearman',...
    'NumDimensions',2) ; % ,'Perplexity', 40); % spearman
% ==================================================================== % 

color = [0 0.4509 0.7411; 0.8705882 0.490196 0];
figure()
gscatter(yValues(:,1),yValues(:,2), groups,color,'..',10)
set(gca,'LineWidth',1,'Box','off')
xlabel('tSNE-1') ; ylabel('tSNE-2') ; legend({'LM','HM'});
title('mRNA-based Grouping of Tumours','FontSize',14','FontWeight','bold')

rng('default') % for reproducibility 
colors = rand(length(unique(groupsCancer)),3);
figure()
gscatter(yValues(:,1),yValues(:,2), groupsCancer,colors,'..',10)
set(gca,'LineWidth',1,'Box','off')
xlabel('tSNE-1') ; ylabel('tSNE-2')
title('mRNA-based Grouping of Tumours','FontSize',14','FontWeight','bold')

% ==================================================================== % 
fprintf('\n Now Running tSNE for 3 Dimesions \n')
% use maybe spectrial clustering here
rng(1) % for reproducibility
yValues3 = tsne(mrnaExpr','Algorithm','barneshut','Distance','spearman',...
    'NumDimensions',3) ; % ,'Perplexity', 40); % spearman
% ==================================================================== % 

% get the colors for the plots to be produced. first get colors for the 3D
% plot for the groups 
color3 = zeros(size(yValues,1),3);
color3(groups,:) = repmat(color(2,:), sum(groups), 1) ;
color3(~groups,:) = repmat( color(1,:), sum(~groups), 1);

% produce the 3D plot 
figure()
scatter3(yValues3(:,1),yValues3(:,2),yValues3(:,3),10,color3,'filled')
set(gca,'LineWidth',1)
xlabel('tSNE-1') ; ylabel('tSNE-2');zlabel('tSNE-3')
title('3D Plot of mRNA-based Grouping of Tumours','FontSize',14',...
    'FontWeight','bold')

% and colors for the 3D
rng(2)
color3Cancer = zeros(size(yValues,1),3);
myCancers = cellstr(unique(groupsCancer)) ;
for ii = 1:length(myCancers)
    color3Cancer(ismember(cellstr(groupsCancer),myCancers(ii)), : ) = ...
        repmat( colors(ii,:), ...
        sum(ismember(cellstr(groupsCancer),myCancers(ii))) ,1);
end

% produce the 3D plot of the cancer types 
figure()
scatter3(yValues3(:,1),yValues3(:,2),yValues3(:,3),10,color3Cancer,'filled')
set(gca,'LineWidth',1)
xlabel('tSNE-1') ; ylabel('tSNE-2');zlabel('tSNE-3')
title('3D Plot of mRNA-based Grouping of Tumours','FontSize',14',...
    'FontWeight','bold')

% ============== create a table and save the data to excel  ==============

tableTSNE = table(yValues(:,1), yValues(:,2), groups, groupsCancer,...
    'VariableNames',{'tSNE1','tSNE2','MetabolicStatus','CancerType'}) ;

writetable(tableTSNE,'t-SNE for Tableau Plotting.xlsx','Sheet','TCGA tSNE');

% also save a table for the gdsc drugs 
writetable(gdscDrugs ,'gdscDrugsTableau.xlsx','Sheet','gdscDrugs');

% check for differentially expressed genes between the two MH and ML
% subtypes and report on what the genes mean for the metabolic pathways. 

colorM  = cbrewer('seq', 'RdPu', 100); % BuPu

% identify the columns with nan values 
these = [false, false, false, any(isnan(mrna4Corr{:,4:end}), 1)] ;
mrna4Corr(:,these) = [] ;
mrnaGstats = grpstats(mrna4Corr,'CancerStudy','mean');

%% =============== Find the Inter Tumour Correlation ===================

% delete the row with LIHC which makes the analysis crap 
mrnaGstats(16,:) = [];
mrnaGstats = sortrows(mrnaGstats,'mean_MetabolicStatus','descend');
mrnaCorr = round(corrcoef(mrnaGstats{:,4:end}'),2)' ;
% [~ ,mrnaCorrPvalue ] = corrcoef(mrnaGstats{:,4:end}') ;

% produce the heatmap
figure(); clf;
axes('position',[0.10, 0.1, 0.8, 0.78]);
heatmap(cellstr(mrnaGstats.CancerStudy),...
    cellstr(mrnaGstats.CancerStudy),mrnaCorr,...
    'Colormap',colorM,'ColorbarVisible','on');

% add altra bars of the cancer metabolic status to the top of the heatmap 
axes('position',[0.10,0.89, 0.8, 0.03]);
heatmap( double(mrnaGstats.mean_MetabolicStatus'+1), ...
    'Colormap',color,'CellLabelColor','none',...
    'FontColor','none','ColorbarVisible','off');

% add altra bars of the of cancer tissue of origin to the top of the heatmap 
axes('position',[0.10,0.93, 0.8, 0.03]);
% change the arrangement in the cancer studies
[~, locX] = ismember(mrnaGstats.CancerStudy, ...
    cellstr(cancerStudies.cancerTypeId)) ;
locX(locX == 0) = [] ;
barData = double( cancerStudies.cancerClass(locX) ) ;
% now get the data of the cancer study
heatmap(double(barData)','Colormap',classesColors.cancerClass,...
    'FontColor','none','ColorbarVisible','off');

% ================= Correlation Mean Analysis =========================
% get the mean and range of the correlation between cancer studies of the
% mRNA levels
% first replace the diagnal with NaN because they involved self comparision
mrnaCorr = mrnaCorr - diag(diag(mrnaCorr)) + diag(NaN(1,length(mrnaCorr)));
locMetab = logical(mrnaGstats.mean_MetabolicStatus);
HMmRNA = mrnaCorr(locMetab,locMetab) ;
LMmRNA = mrnaCorr(~locMetab,~locMetab) ;

% finally get the mean and the range
HMmean = nanmean(nanmean(HMmRNA));
display(HMmean)
HMrange = [nanmin(nanmin(HMmRNA)) ,nanmax(nanmax(HMmRNA)) ] ;
HMpercetile = prctile(HMmRNA,[25 75],'all') ;
display(HMrange)
display(HMpercetile)

LMmean = nanmean(nanmean(LMmRNA));
display(LMmean) 
% repalce the nan with the mean and calculate the range
LMrange = [ nanmin(nanmin(LMmRNA)) ,  nanmax(nanmax(LMmRNA)) ] ;
LMpercetile = prctile(LMmRNA,[25 75],'all') ;
display(LMrange)
display(LMpercetile)

clear pureData mask genes groups mrnaExpr toGo colors myCancers ...
    tableTSNE yValues yValues3 color3Cancer groupsCancer color3 ...
    these toGo them times theUnique theLocs nMetabHigh nMetaboLow ...
    nanIndices mask ii jj j heatData genes mrnaCorr mrnaGstats ...
    mrna4Corr colorM HMrange HMmean LMmean LMrange locMetab expr ...
    hDist LMmRNA locX appendTable curCancer p participants PFSvalues ...
    stats terms ...
    themDrugs trialsResults outcomeResults meanMetab meanNonMetab ...
    locTerm ii kk curDrug curClinTrial dataSpec data ...
    cur_NCTnumber centrals arranGment armLoc ans aa reactomeIDs ...
    DoseResponse colorM 

%% Get the TCGA Classification of Cell Lines and Plot a Graph

% Here the check the metabolic pathways alteratios of cell lines that
% belong to the same cancer and look at the differences in the mutations
% profile and the expresssion of mRNA within the metabolic pathways to see
% if they affect the drug response profile

% first get the number of GDSC cell lines per cancer study
gdscCNA.CancerStudy = categorical(gdscCNA.CancerStudy);
gdscCellLines = grpstats(gdscCNA(:,1),'CancerStudy')  ;
gdscCellLines.Properties.VariableNames(1) = "cancerTypeId" ;
% bar(gdscCellLines.CancerStudy, gdscCellLines.GroupCount)

% produce grouped bar groups using tableau that show the type of cance
% study and TCGA classification and Do the same for the TCGA cancer Studies
% data NO NEED for the sunburst chart.

% and get the misssing cancer descriptions from the gdscCell lines table
% that is saved in the excel file
gdscCellLines = outerjoin(gdscCellLines, cancerStudies(:,[1,2,4]),...
    'MergeKey',true);

% plot box plots of the most responsive drugs based on the pathways that
% the drugs targets: these are EGFRSignaling and ERKMAPKSignaling and
% Metabolism

% save data to excel 
writetable(cancerStudies, 'Supplementary File.xlsx','Sheet',...
    'TCGA Cancer Studies')

writetable(metabolicPathways, 'Supplementary File.xlsx','Sheet',...
    'Metabolic Pathways - First Tier')

% writetable(mutations, 'Supplementary File.xlsx','Sheet',...
%     'Metabolic Pathway Mutations')
% 
% writetable(cnaData, 'Supplementary File.xlsx','Sheet',...
%     'Metabolic Pathway Copy Number')
% 
% writetable(mrna, 'Supplementary File.xlsx','Sheet',...
%     'Metabolic Pathways mRNA')

writetable(gdscCellLines, 'Supplementary File.xlsx','Sheet',...
    'GDSC Cancer Studies')

%% Compare Mutations Between Cell Line and Primary Tumours

% preallocation 
GDSC_TCGA_mutCompare = cellstr(commonCancer) ;
for ii = 1:length(commonCancer)
    % get the rows that have the cancer data
    curRows = ismember( strrep(TCGA_GDSCpathways.CancerStudy,'cl',''),...
        cellstr(commonCancer(ii))  ) ;
    
    % get the total percentage alterations for that gene and add the
    % percentage not altered to the second column of curAlters
    curAlters = round( mean( TCGA_GDSCpathways{curRows,2:end} , 2) );
    [~,p,chi] = prop_test(curAlters(:,1)', [100, 100], true ) ;
    
    % add the values to the comparison table 
    GDSC_TCGA_mutCompare(ii, 2:3) = num2cell([chi,p]) ; 
end
% covert to table 
GDSC_TCGA_mutCompare = array2table(GDSC_TCGA_mutCompare) ;
GDSC_TCGA_mutCompare{:,4} = ...
    mafdr(cell2mat(GDSC_TCGA_mutCompare{:,3}),'BHFDR',true);
GDSC_TCGA_mutCompare.Properties.VariableNames = ...
    {'CancerStudy','ChiSquare','pValue','AdjPValue'} ;
GDSC_TCGA_mutCompare = sortrows(GDSC_TCGA_mutCompare,'AdjPValue','ascend')

writetable(GDSC_TCGA_mutCompare, 'Supplementary File3.xlsx','Sheet',...
    'GDSC-TCGA Alterations')

% clear some variables as always
clear p chi ii curRows curAlters

%% Find Percentage Alteration in each Central Metabolic Pathway

centralTotals = rows2vars(centralPathwaysAlterations,...
    'VariableNamesSource','CancerStudy'  ) ;
centralTotals = addvars(centralTotals(:,1), ...
    sum(centralTotals{:,2:end},2)/(width(centralTotals)-1) );

%% Validate My Find that HM tumours have Worse Survival 

% Here I use ICGC data

% load the data ICGC clinical data and TCGA conversion table
icgc = readtable('icgcValidationSet.txt');
icgctcgaCon = readtable('ICGC to TCGA.xlsx');

% remove TCGA data from the analysis
icgc  = icgc(~contains(icgc.submitted_donor_id, 'TCGA') ,:) ;
icgc.project_code = extractBefore(icgc.project_code ,'-');

% also remove the TCGA 

% return only the the tumours that match with TCGA for both tables to used
% in the comparison
them = ismember(icgctcgaCon.TCGAnames, cancerStudies.cancerTypeId );
icgctcgaCon = icgctcgaCon(them,:);

icgc = icgc( ismember(icgc.project_code, icgctcgaCon.TCGAnames ) ,:);

icgc = addvars(icgc, ...
    ~ismember(icgc.project_code, lowMetabolismTumourNames),...
    'NewVariableNames','Metabolism','Before','project_code');

% perform survival analysis on these independent data

% Get the Overall Survival Data and the Disease Free Survival Data and
% delete the missing low from the data
OsData = [num2cell(icgc.donor_survival_time), ...
    icgc.donor_vital_status ,num2cell(icgc.Metabolism)] ;
OsData(any(cellfun(@isempty,OsData),2),:) = [] ;

fprintf('\nNumber of Low Metabolism Tumours = %d \n',...
    sum(cell2mat(OsData(:,3)) ))
fprintf('\nNumber of High Metabolism Tumours = %d \n',...
     sum(~cell2mat(OsData(:,3)) ))

% Anonotate the groups to either high metabolism or low metabolism and
% perferom K-M Survival Analysis 
groups = cellstr( num2str(cell2mat(OsData(:,3)) ) ) ;
groups = strrep(groups,'0','LM') ;
groups = strrep(groups,'1','HM') ;

% ============= PERFORM OVERALL SURVIVAL ANALYSIS =================
[~, ~, stats] = MatSurv( cell2mat(OsData(:,1)), ...
    lower(OsData(:,2)) , groups, 'NoRiskTable',true) ;

% annotate the plot
set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold',...
    'TickDir','out')
ylevel = 0.65 ;
annotationData = cell2mat(OsData(:,1));
for jj = 1:numel(unique(groups))
    text(max(annotationData)/3,ylevel, ...
        sprintf('Median OS: (%s) = %g\n', ...
        stats.GroupNames{jj},stats.MedianSurvivalTime(jj)),...
        'FontWeight','bold','FontSize',10)
    ylevel = ylevel+0.04 ;
end  
xlabel('Time (Days)')
title('\bf Overall Survival of ICGC Patients','FontSize',16)

% Clear some Variables
clear stats ylevel jj groups icgctcgaCon locId OsData x y yt ytl ...
    barData lineData nMuts nMetabLow nNonMuts them

fprintf('\n THE END!!! \n')

