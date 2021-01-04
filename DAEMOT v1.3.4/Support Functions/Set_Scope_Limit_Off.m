open_bd = find_system('BlockType','Scope');
count = 0;
for i=1:length(open_bd)
    value = get_param(open_bd{i},'LimitDataPoints');
    if strcmp(value,'on')
        count = count+1;
    end
    set_param(open_bd{i},'LimitDataPoints','off')
end
if count == 1
    disp(['Fixed ',num2str(count),' scope!'])
else
    disp(['Fixed ',num2str(count),' scopes!'])
end
