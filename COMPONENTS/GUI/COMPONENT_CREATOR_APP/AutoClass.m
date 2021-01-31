% class heading until properties
clear all
clc

%% Inputs
Name = 'Tank';
Props = {"fluid" "Mass" "Initial_Temp"};
Val = {'H20' '100' '25'};
Com = {'Fluid' 'Mass [kg]' 'Inital Temp [C]'};
CusCom = 'The Tank model stores fluid as a thermal mass';

BuildCompClass(Name,Props,Val,Com,CusCom);




