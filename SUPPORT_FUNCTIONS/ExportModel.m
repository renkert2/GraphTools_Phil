function ExportModel(Name,Model,opt)
disp(['Exporting model ',Name,' to python...'])
path = pwd;
comp = Name;
Model.Nu = height(Model.InputTable);

% create folder
if not(isfolder(comp))
       mkdir(comp)
end

tic
% get jacobians
Ju = jacobian(Model.g_sym,Model.SymVars.u);
Jd = jacobian(Model.g_sym,Model.SymVars.d);

% save tables
saveLoc = [path,'\',comp,'\'];
writetable(Model.DisturbanceTable,[saveLoc,'Dist.csv'])
writetable(Model.InputTable,[saveLoc,'Inp.csv'])
writetable(Model.OutputTable,[saveLoc,'Out.csv'])

% save python files
disp('Working on CalcX ...')
my_mat2py('X','CalcX',Model.g_sym,{Model.SymVars.x Model.SymVars.u Model.SymVars.d},opt)
if Model.Nu > 0
    disp('Working on CalcJu ...')
    my_mat2py('Ju','CalcJu',Ju,{Model.SymVars.x Model.SymVars.u Model.SymVars.d},opt)
end
disp('Working on CalcJd ...')
my_mat2py('Jd','CalcJd',Jd,{Model.SymVars.x Model.SymVars.u Model.SymVars.d},opt)

% my_fprintMatPy2('X',{'x','u','d'},Model.g_sym,0,[length(Model.SymVars.x),length(Model.SymVars.u),length(Model.SymVars.d)])
% my_fprintMatPy2('Ju',{'x','u','d'},Ju,0,[length(Model.SymVars.x),length(Model.SymVars.u),length(Model.SymVars.d)])
% my_fprintMatPy2('Jd',{'x','u','d'},Jd,0,[length(Model.SymVars.x),length(Model.SymVars.u),length(Model.SymVars.d)])

movefile('X.py',comp)
if Model.Nu > 0
    movefile('Ju.py',comp)
end
movefile('Jd.py',comp)
t = toc;
disp(['Export Time: ',num2str(t)]) 
fprintf('\n')

end

