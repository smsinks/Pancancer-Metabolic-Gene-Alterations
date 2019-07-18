% Process mutation data from a MAF into a table 
function mutationTable = processMAF(MAFfile)

% get the CCLE ids and the genes in the table
caseIDs = unique(MAFfile.Tumor_Sample_Barcode) ;
allGenes = unique(MAFfile.Hugo_Symbol)' ;
for ii = 1:length(caseIDs)
    % get all the mutated genes for that sample id
    curMutations = MAFfile( strcmp( ...
        MAFfile.Tumor_Sample_Barcode, caseIDs{ii}),:) ;
    
    % remove duplicate gene mutations as they causing a lot of error and
    % also remove the silent mutations from the table
    curMutations(strcmp(curMutations.Variant_Type,...
        'Silent'),:) = [] ;
    [~,ia] = unique(curMutations(:,1),'rows','first');
    curMutations = curMutations(ia,:);

    % now compare the curMutations with allGenes this comparison is only
    % for the first row
    [~, posLocs]= intersect(allGenes(1,:),curMutations.Hugo_Symbol ) ;
    
    % all to the genes list 
    allGenes(ii+1,posLocs') = curMutations.Variant_Type ; 
    
end

% now combine the two cell array and convert in to a table
mutationTable = cell2table([caseIDs, allGenes(2:end,:)] ) ;
mutationTable.Properties.VariableNames = ...
    matlab.lang.makeValidName ( ['SampleIDs',allGenes(1,:)]  );

% delete all the genes that do not have mutations in the table
toGo = all( cellfun(@isempty,mutationTable{:,:} ) );
mutationTable(:,toGo) = [] ;