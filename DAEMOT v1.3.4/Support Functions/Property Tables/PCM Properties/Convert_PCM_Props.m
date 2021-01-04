%   Save PCM Properties into PCMProp Structure
load('PCMProps.mat')
names = fieldnames(PCMProps);

[~,~,A]=xlsread('PCM Thermophysical Properties.xlsx');

k=0;
for n=4:length(A(2,:))
    for n2=1:length(names)
        k=k+strcmp(A(2,n),names{n2});
    end
    
    if k==0%new material name has been entered
        new_ID = A(2,n);
        new_ID=new_ID{1};
        new_props=[A(3:9,1) A(3:9,n)];
        new_props = cell2struct(new_props(:,2),new_props(:,1),1);
        PCMProps.(new_ID)=new_props;
    end
    k=0;
end

%Uncomment the following if saving new values
%save('PCMProps.mat','PCMProps')

%Edit PCM Simulink block mask to update these changes:
   % in the Parameters & Dialog tab, edit the Phase Change Material popup:
   % Click on the Phase Change material popup
   % Edit 'Type options'(popup options)
   % Add the MaterialName specified in the Excel spreadsheet to the list of
   % options

    % in the Initialization tab, add the following to the first if
    % statement (change MaterialName to the name in Excel spreadsheet):
    %elseif strcmp(PCM_ID,'MaterialName')
    %PCM_info=PCMProps.MaterialName;