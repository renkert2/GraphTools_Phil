function obj = my_inputparser(obj,varargin)
    % NameValueDef allows the population of a class or structure using
    % name-value pairs. 
    % 
    
    %%% INPUTS
    % obj - structure/object to be populated
    % varargin - name-value pair information

    %%% OUTPUTS
    % obj - updated structure/object
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Christopher T. Aksland
    % Association: University of Illionis at Urbana-Champaign
    % Contact: aksland2@illinois.edu
    % Revision History:
    % 11/22/2020 - Class creation (Referenced YALMIP sdpsettings function
    %             for development).
    % 11/24/2020 - matlab has an input parser function/class that could
    %             replace this custom function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NAMES = fieldnames(obj); % get field names for the object

if mod(nargin-1,2) ~= 0 % error check to make sure there are an even number of arguements
    error('Provide arguments in name-value pairs.');
end
i = 1; % start indexing at location 1
val = 0; % 0: indicates argument name, 1: indicates argument value

while i <= nargin-1 % loop through all arguements
    par = varargin{i}; % extract the argument
    
    if val == 0 % check if the arguement is a name or value
        if ~ischar(par) % error check if argument is not a string name
            error('Provide name arguement as a string name.')
        end
        
        if ~(startsWith(par,NAMES)) % check if the provided argument is valid
            error([sprintf('Provide a valid argument name:\n') sprintf('-%s\n',NAMES{:})]) % update this to list valid arguments
        end
        val = 1; % denotes that the next argument is a value
        
    else
        obj.(varargin{i-1}) = par; % set the argument value
        val = 0; % denotes that the next argument is a name
        
    end
    i = i + 1; % index through the loop
    
end


end