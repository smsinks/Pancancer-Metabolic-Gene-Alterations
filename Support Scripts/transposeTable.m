% Transpose any table with this function 
%
% Input = any table
%
% Output = a transpose of the input table

function mytable = transposeTable(in_table)

myArray = table2cell(in_table(:,2:end) );
myArray = cell2table(myArray'); 
var_names = cellstr( table2cell(in_table(:,1)) );
var_names = matlab.lang.makeValidName(var_names) ;
var_names = var_names';

myArray.Properties.VariableNames = var_names ;

% S = {'my.Name','my_Name','my_Name'};
% validValues = matlab.lang.makeValidName(S)
% validUniqueValues = matlab.lang.makeUniqueStrings(validValues,{},...
%     namelengthmax)
  
row_names = in_table.Properties.VariableNames(2:end); 
row_names = cell2table(row_names');
mytable = [row_names , myArray ] ;
mytable.Properties.VariableNames(1,1) = in_table.Properties.VariableNames(1,1);

clear myArray var_names row_names ii str expression replace newStr 
end