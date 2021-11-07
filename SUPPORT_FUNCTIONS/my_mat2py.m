function my_mat2py(filename,functname,symFunc,inputs,opt)

%% generate matlab function
symCell = reshape(sym2cell(symFunc),[],1);
N = numel(symCell);
oName = "out"+string(1:N);
x = inputs{1};
u = inputs{2};
d = inputs{3};
matlabFunction(symCell{:},'Outputs',cellstr(oName),'Vars',{x,u,d},'File',functname,'Optimize',opt,'Comments','None');

%%
CalcText = textscan(fopen([functname,'.m']),'%s','Delimiter','\n');
CalcText = string(CalcText{1});

% Replace Header
CalcText(1:8) = [];
funcDef = ["import math"];
funcDef(end+1,1) = [""];
funcDef(end+1,1) = ["def " + functname + "(x, u, d):"];
funcDef(end+1,1) = "# auto-generated function from matlab";
funcDef(end+1,1) = [""];

% replace input variable parsing
CalcText(contains(CalcText,'in')) = [];
dVars = [string(d)+"=d["+string((0:(length(d)-1))')+"]"];
uVars = [string(u)+"=u["+string((0:(length(u)-1))')+"]"];
ParseText = [dVars;uVars;""];

% replace matlab syntax
CalcText = strrep(CalcText,'.*','*');
CalcText = strrep(CalcText,'./','/');
CalcText = strrep(CalcText,'.^','**');
CalcText = strrep(CalcText,'sqrt','math.sqrt');
CalcText(contains(CalcText,'nargout')) = [];
CalcText(contains(CalcText,'end')) = [];

% function return string
returnFunc = ['return ', sprintf('%s, ',oName{1:end-1}),oName{end}];
CalcText(end+1) = returnFunc;

% write to file 
fid = fopen([filename,'.py'],'w');
fprintf(fid,'%s\n',funcDef(:));
fprintf(fid,'\t%s\n',ParseText(:));
fprintf(fid,'\t%s\n',CalcText(:));
fclose(fid);

end