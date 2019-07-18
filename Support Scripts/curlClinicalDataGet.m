%% curlClinicalDataGet
% Get clinical data from cBioPortal using a Curl API
% get the location of the data on the cBioPortal website
% dataLocation = sprintf(['http://www.cbioportal.org/api/studies/%s/',...
%     'clinical-data?clinicalDataType=PATIENT&projection=',...
%     'SUMMARY&pageSize=10000000&pageNumber=0&direction=ASC'],cancerTypeId);
%    tempTable = struct2table( webread(dataLocation));

cancerStudies = readtable('cancerStudies.txt');
cancerStudies.Properties.VariableNames(1) = {'cancerTypeId'};

% Get the list cancer studies codes & return only Pancaner TCGA studies for
% whcih the data has been normalise for ease of comparison
cancerStudies = cancerStudies( contains(cancerStudies.name, ...
    'PanCancer','IgnoreCase',true) ,:) ;

%% Access the data without Using Curl

clinicalData = [] ;
for ii = 1:height(cancerStudies)
    % get the current cancer study
    cancerTypeId = string ( extractBefore( ...
        cancerStudies.cancerTypeId(ii), '_pan')) ;
    
    fprintf('Now Processing Clinical Data for: %s Study # %d\n',...
        cancerTypeId,ii)
    
    getDataName = strcat(cancerTypeId,'_sequenced') ;
    tempTable = getclinicaldata(cgdsURL, sprintf('%s',getDataName));
    appendTable = array2table([tempTable.caseId, tempTable.data]) ;
    appendTable.Properties.VariableNames = ['SampleIds';...
        tempTable.clinVariable ] ;
%     % convert the current table into a proper table
%     patientIds = unique(tempTable.patientId) ;
%     appendTable = unique(tempTable(:,4)) ;
%     for jj = 1:length(patientIds)
%         % get the data for the current patient
%         patientLoc = tempTable( ismember(tempTable.patientId,...
%             patientIds(jj)),[4,5]);
%         patientLoc.Properties.VariableNames(2) = ...
%             strrep(patientIds(jj),'-','_') ;
%         
%         % append the patient to the growing table
%         appendTable = innerjoin(appendTable,patientLoc);
%     end
    
   
%     appendTable = rows2vars(appendTable,'VariableNamesSource',...
%         'clinicalAttributeId');
%     appendTable.Properties.VariableNames(1) = {'SampleIds'};
    
    % add cancer study
    appendTable = addvars(appendTable,...
        repmat( cancerTypeId, height(appendTable),1) ,...
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
            % [C,ia,ib] = intersect(___) also returns index vectors ia
            % and ib using any of the previous syntaxes. Generally, C =
            % A(ia) and C = B(ib)
            %             [~,ia,ib] = intersect(clinicalData.Properties.VariableNames,...
            %                 appendTable.Properties.VariableNames,'stable') ;
            %             clinicalData = clinicalData(:,ia);
            %             appendTable = appendTable(:,ib);
            
            % sometimes the table variables are cell array and not double
            % therefore I have to take that into account
            clinicalData = vertcat(clinicalData, appendTable) ;
    end
end