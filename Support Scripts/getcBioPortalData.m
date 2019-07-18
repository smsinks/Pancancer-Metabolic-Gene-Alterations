
% This function processes TCGA data for pancancer studies into single files
% for each dataset for all studies

function [mutations, cancerStudies, cnaData, mrna, clinicalData ] = ...
    getcBioPortalData(myGenes)

% the api seem to be non functional at the moment. So I use a text filed
% that i download from cBioportal to run my analysis
try
    % Get Data from cBioPortal
    % Set web API URL (excluding 'webservice.do', trailing slash optional)
    cgdsURL = 'http://www.cbioportal.org/';
    % Get list of available cancer types
    cancerStudies = getcancerstudies(cgdsURL);
    cancerStudies = struct2table(cancerStudies);
catch
    cancerStudies = readtable('cancerStudies.txt');
    cancerStudies.Properties.VariableNames(1) = {'cancerTypeId'};
end

% Get the list cancer studies codes & return only Pancaner TCGA studies for
% whcih the data has been normalise for ease of comparison
cancerStudies = cancerStudies( contains(cancerStudies.name, ...
    'PanCancer','IgnoreCase',true) ,:) ;


%% Process the Copy Number Data into One Table

% This fetches mutation, copy number alterations and expression data for
% all the genes (allGenes) for each patient ID

cnaData = [] ; % this also helps to clear the present CNAdata
for ii = 1:height(cancerStudies)
    % read the cancer study from the text files one after the other. The
    % cancer Ids are found in the cBioPortal table "cancerStudies"
    cancerTypeId = string ( extractBefore( ...
        cancerStudies.cancerTypeId(ii), '_tcga')) ;
    tempCNA = readtable( strcat(cancerTypeId,'_tcga_data_CNA.txt') );
    fprintf('Now Processing Copy NumberData for: %s Study # %d\n',...
        cancerTypeId,ii)
    
    % get the genes which are found in metabolic pathways. Also remove the
    % Entrez Gene Ids and then transponse the table so that the genes are
    % not top
    tempCNA = tempCNA( ismember(tempCNA.Hugo_Symbol,myGenes) ,: ) ;
    tempCNA.Entrez_Gene_Id = [] ;
    
    % also some of the data have cytobands of chromosome. delete these
    try tempCNA.Cytoband = [] ; catch ; end
    
    tempCNA = rows2vars(tempCNA,'VariableNamesSource','Hugo_Symbol');
    tempCNA.Properties.VariableNames(1) = {'SampleIds'};
    tempCNA = addvars(tempCNA,...
        upper(repmat(cancerTypeId,height(tempCNA),1)) ,...
        'Before','SampleIds','NewVariableNames','CancerStudy') ;
    
    % find the common genes between the two dataset and then combine the
    % table.
    switch ii
        case ii == 1
            cnaData = vertcat(cnaData,tempCNA) ;
        otherwise
            %     [C,ia,ib] = intersect(___) also returns index vectors ia
            %     and ib using any of the previous syntaxes. Generally, C =
            %     A(ia) and C = B(ib)
            [~,ia,ib] = intersect(cnaData.Properties.VariableNames, ...
                tempCNA.Properties.VariableNames,'stable') ;
            cnaData = cnaData(:,ia);
            tempCNA = tempCNA(:,ib);
            
            % sometimes the table variables are cell array and not double
            % therefore I have to take that into account
            try
                cnaData = vertcat(cnaData, tempCNA) ;
            catch
                fprintf('The Data %d is of Cell Type\n',ii)
                for jj = 3:width(tempCNA)
                    tempCNA.(jj) = str2double(tempCNA.(jj)) ;
                end
                cnaData = vertcat(cnaData, tempCNA) ;
            end
    end
end

% remove the 1 and - 1 or - 0
for ii = 3:width(cnaData)
    cnaData.(ii)(cnaData.(ii) == -1|cnaData.(ii) == 1|cnaData.(ii) == -0) = 0  ;
end

clear tempCNA jj cancerTypeId ia ib

%% Process The Clinical and Sample Data into a Single Table

clinicalData = [] ;
fprintf('\n\n')
for ii = 1:height(cancerStudies)
    % get the current cancer study
    cancerTypeId = string ( extractBefore( ...
        cancerStudies.cancerTypeId(ii), '_pan')) ;
    
    fprintf('\nNow Processing Clinical Data for: %s Study # %d\n',...
        cancerTypeId,ii)
    
    getDataName = strcat(cancerTypeId,'_sequenced') ;
    tempTable = getclinicaldata(cgdsURL, sprintf('%s',getDataName));
    appendTable = array2table([tempTable.caseId, tempTable.data]) ;
    appendTable.Properties.VariableNames = ['SampleIds';...
        tempTable.clinVariable ] ;
    
    % add cancer study
    appendTable = addvars(appendTable,...
        upper( repmat( cancerTypeId, height(appendTable),1) ),...
        'Before','SampleIds','NewVariableNames','CancerStudy') ;
    
    % this is different from the copy number data as not all the
    % clinical data . Therefore I need to add the genes that are
    % have the same names the table create a table of
    if ii ~=1
        missingClinicals = setdiff(...
            clinicalData.Properties.VariableNames, ...
            appendTable.Properties.VariableNames) ;
        dummyClins1 = array2table( cell(height(appendTable),...
            length(missingClinicals)) );
        dummyClins1.Properties.VariableNames = missingClinicals ;
        appendTable = [appendTable,dummyClins1] ;
        
        % also add the dummy clinical data to the growing table
        missingClinicals2 = setdiff(appendTable.Properties.VariableNames,...
            clinicalData.Properties.VariableNames) ;
        dummyClins2 = array2table( cell(height(clinicalData),...
            length(missingClinicals2)) );
        dummyClins2.Properties.VariableNames = missingClinicals2 ;
        clinicalData = [clinicalData,dummyClins2] ;
    end
    
    % add the current clinical data to the overall clinical data
    switch ii
        case ii == 1
            clinicalData = vertcat(clinicalData,appendTable) ;
        otherwise           
            % sometimes the table variables are cell array and not double
            % therefore I have to take that into account
            clinicalData = vertcat(clinicalData, appendTable) ;
    end
end

%% Process the Mutations Data into a Single Table

mutations = [] ; % this also helps to clear the present Mutations data
cancerMuts = zeros(height(cancerStudies),2);
fprintf('\n\n')
for ii = 1:height(cancerStudies)
    % read the cancer study from the text files one after the other. The
    % cancer Ids are found in the cBioPortal table "cancerStudies"
    cancerTypeId = string ( extractBefore( ...
        cancerStudies.cancerTypeId(ii), '_tcga')) ;
    tempMut = readtable( strcat(cancerTypeId,...
        '_tcga_data_mutations_extended.txt') );
    fprintf('Now Processing Muations Data for: %s Study # %d\n',...
        cancerTypeId,ii)
    
    % get teh total number and mutations and them to the table
    fprintf('  The total number of genes is %d \n',...
        length( unique(tempMut.Hugo_Symbol) ) )
    cancerMuts(ii,1) = length( unique(tempMut.Hugo_Symbol) );
    
    % get the genes which are found in metabolic pathways and add them to
    % the cell array 
    tempMut = tempMut( ismember(tempMut.Hugo_Symbol, myGenes) ,: ) ;
    fprintf('  Remaining (Metabolic Pathway) genes is %d\n\n',...
        length( unique(tempMut.Hugo_Symbol) ) )
    cancerMuts(ii,2) = length( unique(tempMut.Hugo_Symbol));
    
    % process the resultant MAF file
    tempMut = processMAF(tempMut);
    
    % change the name from SampleIDs to SampleIds to make the name uniform
    % with other data. Also replace - with _
    tempMut.Properties.VariableNames(1) = {'SampleIds'};
    tempMut.SampleIds = strrep(tempMut.SampleIds,'-','_') ;
    
    % add cancer study
    tempMut = addvars(tempMut, ...
        upper( repmat(cancerTypeId,height(tempMut),1)) ,...
        'Before','SampleIds','NewVariableNames','CancerStudy') ;
    
    % this is different from the copy number data as not all the metabolic
    % pathway genes are mutated. Therefore I need to add the genes that are
    % not mutated to the table create a table of height = that of tempMut
    % and length = to that of the unmutated genes
    missingGenes = setdiff(myGenes,tempMut.Properties.VariableNames(3:end));
    dummyMut = array2table( cell(height(tempMut),length(missingGenes)) );
    dummyMut.Properties.VariableNames = strrep(missingGenes,'-','_') ;
    tempMut = [tempMut,dummyMut] ;
    
    % add to the growing table
    switch ii
        case ii == 1
            mutations = vertcat(mutations,tempMut) ;
        otherwise
            %     [C,ia,ib] = intersect(___) also returns index vectors ia
            %     and ib using any of the previous syntaxes. Generally, C =
            %     A(ia) and C = B(ib)
            [~,ia,ib] = intersect(mutations.Properties.VariableNames, ...
                tempMut.Properties.VariableNames,'stable') ;
            mutations = mutations(:,ia);
            tempMut = tempMut(:,ib);
            mutations = vertcat(mutations, tempMut) ;
    end
end

cancerStudies = addvars(cancerStudies, cancerMuts(:,1), cancerMuts(:,2),...
    'NewVariableNames',{'totalMutations','metabolicMutations'}) ;

% remove the duplicate rows of samples Ids
[~, nonDuplicates] = unique(mutations.SampleIds);
mutations = mutations(nonDuplicates,:);

% clear some variables
clear tempMut ia ib ii cancerTypeId dummyMut missingGenes

%% Get the Gene Expression Profile of the Metabolic Genes

mrna = [] ; % this also helps to clear the present Mutations data
for ii = 1:height(cancerStudies)
    % read the cancer study from the text files one after the other. The
    % cancer Ids are found in the cBioPortal table "cancerStudies"
    cancerTypeId = string ( extractBefore( ...
        cancerStudies.cancerTypeId(ii), '_tcga')) ;
    tempRNA = readtable( strcat(cancerTypeId,...
        '_tcga_data_RNA_Seq_v2_expression_median') );
    fprintf('Now Processing mRNA Data for: %s Study # %d\n',cancerTypeId,ii)
    fprintf('  The total number of genes is %d \n',...
        length( unique(tempRNA.Hugo_Symbol) ) )
    
    % get the genes which are found in metabolic pathways.
    tempRNA = tempRNA( ismember(tempRNA.Hugo_Symbol, myGenes) ,: ) ;
    fprintf('  Remaining genes is %d\n\n',...
        length( unique(tempRNA.Hugo_Symbol) ) )
    % remove the entrenz gene symbol and convert the cell array data in the
    % table to double
    tempRNA.Entrez_Gene_Id = [] ;
    try
        tempRNA{:,3:end} = strrep(tempRNA{:,2:end},'NA','NaN');
    catch
    end
    % change the cell containts to double
    for jj = 2:width(tempRNA)
        if iscell(tempRNA.(jj))
            tempRNA.(jj) = str2double(tempRNA.(jj)) ;
        end
    end
    % transpose the table
    tempRNA = rows2vars(tempRNA,'VariableNamesSource','Hugo_Symbol') ;
    tempRNA.Properties.VariableNames(1) = {'SampleIds'};
    
    % add cancer study
    tempRNA = addvars(tempRNA,repmat(cancerTypeId,height(tempRNA),1) ,...
        'Before','SampleIds','NewVariableNames','CancerStudy') ;
    
    % this is different from the copy number data as not all the metabolic
    % pathway genes are mutated. Therefore I need to add the genes that are
    % not mutated to the table create a table of height = that of tempMut
    % and length = to that of the unmutated genes
    missingGenes =setdiff(myGenes,tempRNA.Properties.VariableNames(3:end));
    dummyRNA = array2table( cell(height(tempRNA),length(missingGenes)) );
    dummyRNA.Properties.VariableNames = strrep(missingGenes,'-','_') ;
    tempRNA = [tempRNA,dummyRNA] ;
    
    % add to the growing table
    switch ii
        case ii == 1
            mrna = vertcat(mrna,tempRNA) ;
        otherwise
            %     [C,ia,ib] = intersect(___) also returns index vectors ia
            %     and ib using any of the previous syntaxes. Generally, C =
            %     A(ia) and C = B(ib)
            [~,ia,ib] = intersect(mrna.Properties.VariableNames, ...
                tempRNA.Properties.VariableNames,'stable') ;
            mrna = mrna(:,ia);
            tempRNA = tempRNA(:,ib);
            mrna = vertcat(mrna, tempRNA) ;
    end
end

% clear some variables
clear tempMut ia ib ii jj cancerTypeId dummyRNA missingGenes

%% Get the common samples for each all genetic profiles

% extract the sample Ids before the -01. This is what is making the
% difference. First replace all end of the name with _01
mrna.SampleIds = strcat( extractBefore(mrna.SampleIds,13), '_01');
cnaData.SampleIds = strcat( extractBefore(cnaData.SampleIds,13), '_01');
mutations.SampleIds = strcat( extractBefore(mutations.SampleIds,13),'_01');
clinicalData.SampleIds = strrep( strcat( extractBefore(...
    clinicalData.SampleIds,13),'_01') , '-', '_');

commonIds = intersect(cnaData.SampleIds,mutations.SampleIds);
fprintf('\nThe number of Pancancer Patients is %d \n',length(commonIds) )

% retain only the row with matching IDs
mrna = mrna( ismember(mrna.SampleIds, commonIds), :) ;
cnaData = cnaData( ismember(cnaData.SampleIds, commonIds), :) ;
mutations = mutations( ismember(mutations.SampleIds, commonIds), :);
clinicalData = clinicalData(ismember(clinicalData.SampleIds,...
    commonIds),:);

% convert the cancerType to categories
cancerStudies.cancerTypeId = ...
    categorical(upper(extractBefore(cancerStudies.cancerTypeId ,'_tcga')));


end % end of function
