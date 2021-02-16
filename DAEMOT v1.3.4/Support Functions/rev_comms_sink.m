function rev_comms_sink()

h      = get_param(gcb,'handle');   % Get block handle                           
hblock = get(h);     % Parse out handles structure         
str = [hblock.Path, hblock.Name];
str(~ismember(str,['A':'Z' 'a':'z' '0':'9' '_'])) = '';                      
set_param([gcb, '/Goto'],'GotoTag',['gdata_',str]);
                                                                                 
clear h hblock str           