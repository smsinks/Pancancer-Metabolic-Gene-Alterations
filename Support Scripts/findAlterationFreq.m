% This function returns the alteration frequency for mutations and copy
% number data

function pathwayAlterations = findAlterationFreq(metabolicPathways,...
    mutations,cnaData)
% This is a table with columns as cancer types and row as metabolic
% pathways

% get the genes involved in a metabolic pathwy and return only these for
% the copy number data and mutations data

for ii = 1:height(metabolicPathways)
    
    % get the genes involved
    pathwayGenes = split(metabolicPathways.Genes(ii));
    pathwayMuts = mutations(:,[true, false, ...
        ismember(mutations.Properties.VariableNames(3:end), pathwayGenes)]);
    
    % process the copy number add
    if nargin == 3
        pathwayCopyNumber = cnaData(:,[true, false, ...
            ismember(cnaData.Properties.VariableNames(3:end),...
            pathwayGenes)] );
    end
    
    % get the mutatations in each samples
    pathwayMuts = addvars( pathwayMuts(:,1), ...
        double(any(~cellfun(@isempty,pathwayMuts{:,2:end}),2) ) , ...
        'NewVariableNames','Overall') ;
    pathwayMuts.CancerStudy = categorical(pathwayMuts.CancerStudy);
    
    % combine the two tables: ONLY IF there are also copy number
    % alterations
    if nargin == 3
        if any(ismember(cnaData.Properties.VariableNames(3:end),...
                pathwayGenes))
            
            % also get alterations for the copy number data. Sometimes the
            % data is cell if iam deleteion with GDSC data
            try % for TCGA and CCLE data
                pathwayCopyNumber2 = addvars( pathwayCopyNumber(:,1), ...
                    double( any(pathwayCopyNumber{:,2:end}, 2) ) , ...
                    'NewVariableNames','Overall') ;
            catch % for GDSC data
                pathwayCopyNumber2 = addvars( pathwayCopyNumber(:,1), ...
                    double(any(...
                    ~cellfun(@isempty,pathwayCopyNumber{:,2:end}),2) ),...
                    'NewVariableNames','Overall') ;
                
            end
            pathwayCopyNumber2.CancerStudy = ...
                categorical(pathwayCopyNumber2.CancerStudy);
            
            % add to the total mutations table
            pathwayMuts.Overall =  double( ...
                any([pathwayMuts.Overall,pathwayCopyNumber2.Overall] ,2) ) ;
        end
    end
    % convert the zeroes to NaN to make group stats easiler to do
    pathwayMuts.Overall(pathwayMuts.Overall == 0) = NaN ;
    
    % now get the group stats for the two tables: first get the stats for
    % the table without zeros and then with zeroes
    pathwayMuts = grpstats(pathwayMuts,'CancerStudy','numel') ;
%     pathwayCopyNumber =
%     grpstats(pathwayCopyNumber,'CancerStudy','numel');
    
    % create a table that has the over alterations percentage for both
    % mutations and copy number data
    if ii == 1
        pathwayAlterations = addvars(pathwayMuts(:,1), ...
            round( ...
            pathwayMuts.numel_Overall./pathwayMuts.GroupCount,3)*100,...
            'NewVariableNames', ...
            matlab.lang.makeValidName(metabolicPathways.pathwayName(ii)) );
    else % % join the two tables
        tempAlterations = addvars(pathwayMuts(:,1), ...
            round( ...
            pathwayMuts.numel_Overall./pathwayMuts.GroupCount,3)*100,...
            'NewVariableNames', ...
            matlab.lang.makeValidName(metabolicPathways.pathwayName(ii)) );
        
        pathwayAlterations = innerjoin(pathwayAlterations, tempAlterations);
    end
end

end