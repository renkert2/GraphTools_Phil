function Data = TDMS_to_Data_Struct(file,path)
clear Values T Data
file_name=fullfile(path,file);

[ConvertedData]=simpleConvertTDMS(0,file_name);

A = ConvertedData.Data.MeasuredData;

N_channels = size(A,2)-2;

Names = cell(N_channels,1);
for i = 1:N_channels
    T(i) = A(i+2).Total_Samples;
end
Num_Samples = min(T);
for i = 1:N_channels
    Name = A(i+2).Name;
    Name = strsplit(Name,'/');
    Name = strrep(Name, ' ', '_');
    Name = strrep(Name, '(', '');
    Name = strrep(Name, ')', '');
    Names(i,1) = Name(2);
    Values(i,:) = A(i+2).Data(1:Num_Samples);
    Data.(Names{i,1}) = Values(i,:)';
end

