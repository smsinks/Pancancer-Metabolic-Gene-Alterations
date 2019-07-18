% This function returns the alteration frequency for mutations and copy
% number data

function pathwayAlterations = findGeneAlterationFreq(metabolicPathways,...
   studyMuts,studyCNA)

% This is a table with columns as cancer types and row as metabolic
% pathways

% get the genes involved in a metabolic pathway and return only these for
% the copy number data and mutations data

for ii = 1:height(metabolicPathways)
    
    % get the genes involved
    pathwayGenes = split(metabolicPathways.Genes(ii));
    pathwayMuts = studyMuts(:,[false, true, ...
        ismember(studyMuts.Properties.VariableNames(3:end), pathwayGenes)]);
    pathwayCopyNumber = studyCNA(:,[false, true, ...
        ismember(studyCNA.Properties.VariableNames(3:end), pathwayGenes)] );
    
    % get the divisor value to be used to calculate the ratio of pathway
    % alteration in each sample
    divisor = width(pathwayMuts)+ width(pathwayCopyNumber) -2 ;
    
    % get the mutatations in each samples
    pathwayMuts = addvars( pathwayMuts(:,1), ...
        sum(~cellfun(@isempty,pathwayMuts{:,2:end}),2) , ...
        'NewVariableNames','Overall') ;
    
    % also get alterations for the copy number data
    pathwayCopyNumber = addvars( pathwayCopyNumber(:,1), ...
        sum(pathwayCopyNumber{:,2:end} ~= 0, 2) , ...
        'NewVariableNames','Overall') ;
    
    % combine the two tables
    pathwayMuts.Overall = ...
        (pathwayMuts.Overall+pathwayCopyNumber.Overall)/...
        divisor;
    
%     max(pathwayMuts.Overall)
    
    % create a table that has the over alterations percentage for both
    % mutations and copy number data
    if ii == 1
        % add the alteration percentage to the growing table
        pathwayAlterations = addvars(pathwayMuts(:,1), ...
            round( ...
            pathwayMuts.Overall,3)*100,...
            'NewVariableNames', ...
            matlab.lang.makeValidName(metabolicPathways.pathwayName(ii)) );
    else % % join the two tables
        tempAlterations = addvars(pathwayMuts(:,1), ...
            round( ...
            pathwayMuts.Overall,3)*100,...
            'NewVariableNames', ...
            matlab.lang.makeValidName(metabolicPathways.pathwayName(ii)) );
        
        pathwayAlterations = innerjoin(pathwayAlterations, tempAlterations);
    end
end

end