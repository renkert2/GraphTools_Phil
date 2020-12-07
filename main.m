clear all
close all
clc

tank1 = Tank('Name','Tank 1','T_init',50);
load1 = HeatLoad('Name','Load 1','T_init',16);
HX1   = HeatExchanger('Name','HX 1','T1_init',16);
SJ1   = SplitJunction('Name','J3S2','n_in',3,'n_out',2);

Comps = [tank1; load1; HX1; SJ1];


figure
for i = 1:length(Comps)
    subplot(2,2,i)
    plot(Comps(i).Model)
    title(Comps(i).Name)
end
