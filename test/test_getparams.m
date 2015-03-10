%% testing script

clc; clear;

%% testing get_parameters function

% test for wrong number of args
%param_cell0 = get_parameters('s1', 's2');
%param_cell1 = get_parameters('s1', 's2', 's4', 's5');
param_cell2 = get_parameters('var1', '5', 'var2', 'val2', 'var3', '0.5',...
    'var4', 'dd/mm/yy', 'var5', 'NaN', 'user_params.txt');

fprintf('get_parameters test done\n');