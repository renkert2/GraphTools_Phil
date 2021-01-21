function list = generateGraphList()

% this is a function that will need to be updated whenever we add a new
% component the tool functionality. The list name should be the name of the
% class.

i = 1;
list{i} = 'Tank';              i = i + 1;
list{i} = 'HeatExchanger';     i = i + 1;
list{i} = 'HeatLoad';          i = i + 1;
list{i} = 'SplitJunction';     i = i + 1;

list{i} = 'Source';            i = i + 1;
list{i} = 'Sink';              i = i + 1;

end