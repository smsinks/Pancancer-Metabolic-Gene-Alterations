% get the all genes in are in the 3 pathways from the ucsc pathway

function mTOR = createNetwork2(mTOR)


fprintf('\n Reading Super Pathway Data \n')
ucscPathway = readtable('My_Super_Pathway.xlsx');

ogmTOR = mTOR;
mTOR = ucscPathway( contains(ucscPathway.Protein1,mTOR) & ...
    contains(ucscPathway.Protein2, mTOR) , :) ;

% ogCellcycle = cellCycle ;
% cellCycle = ucscPathway( contains(ucscPathway.Protein1,cellCycle) & ...
%     contains(ucscPathway.Protein2,cellCycle) , :);
%
% ogPathwaysInCancer = pathwaysInCancer;
% pathwaysInCancer = ucscPathway( contains(ucscPathway.Protein1,...
%     pathwaysInCancer) & contains(ucscPathway.Protein2,pathwaysInCancer),:);

% ========== These results look very good so far!!!!!!! ============
% clear locK2 locK pathways GOKegg1 GOKegg2 GoKeggC1Terms GoKeggC2Terms ...
%     GoBioC1Terms GoBioC2Terms GObio1 GObio2 to_go combinedScores ...
%     goTerm termEnd ucscPathway

%% Now create the mTOR pathaway

% I will have to start with a simpler (smaller) graph
% create a graph from edge: first delete the selfloop nodes
selfLoop = strcmpi(mTOR.Protein1 , mTOR.Protein2) ;
mTOR(selfLoop,:) = [] ;

% create table for yED that also contains the origanal proteins and also
% add the differential gene expression log p-value to be used to colour the
% notes
ogProteins = contains(mTOR.Protein1 ,ogmTOR) ;
mTOR.NodeInGO = double(ogProteins) ;

% remove that bad arrow from the data
bad = contains(mTOR.Interaction,'-t>');
mTOR.Interaction(bad,1) = {'->t'};

% now create a graph
mtorGraph = digraph(mTOR.Protein1 , mTOR.Protein2);
mtorGraph.Edges.Interaction = mTOR.Interaction ;

% plot the graph
figure()
hMTOR = plot(mtorGraph,'layout','force','usegravity',true,...
    'MarkerSize',10,'ArrowSize', 6,'EdgeAlpha',0.80 ,...
    'LineWidth', 0.5000);
set(gca,'FontSize',12,'FontWeight','bold','visible', 'off')
title('mTOR-Network','FontSize',18)
hold on
% get the nodes that have edge for interactions from biogrid
allInters =  unique(mtorGraph.Edges.Interaction) ;
for ii = 1:length(allInters)
    cur_inter = allInters(ii,1) ;
    locsG = contains(mtorGraph.Edges.Interaction,cur_inter);
    [sOut,tOut] = findedge(mtorGraph);
    allEdges = [sOut,tOut];
    % check = mtorGraph.Edges(locsG,:) ;
    subGraph = allEdges(locsG,:) ;
    subGraph = reshape( subGraph',1,[]) ;
    % if the interaction is just protein-protein
    if strcmp(cur_inter,'->i')
        highlight(hMTOR,subGraph,'EdgeColor',[0.73 0.49 0.43], ...
            'LineWidth',1.5 ,'LineStyle','--') % ,'ArrowPosition',1)
    elseif strcmp(cur_inter,'->p')
        highlight(hMTOR,subGraph,'EdgeColor','b','LineWidth',2)
    elseif strcmp(cur_inter,'-a>')
        highlight(hMTOR,subGraph,'EdgeColor',[0.32 0.79 0.43],...
            'LineWidth',2)
    elseif strcmp(cur_inter,'-a|')
        highlight(hMTOR,subGraph,'EdgeColor','r','LineWidth',2)
    else
        highlight(hMTOR,subGraph,'EdgeColor',[0.5 0.5 0.5],...
            'LineWidth',1.5)
    end
    
end
hold off