function [header] = findspectra(filename)

% Reads output files from starlight and finds the spectrum fit data
% input; filename starlight output file
% output; header = Line in the file at which the fit starts
%         nolines = number of lines that fit takes. Can/should be drop??
%
%     Zachary Byrne


outputfile = fopen(filename, 'r');

%nolines = -1;
outline = fgets(outputfile);
header = 0;
oline = 1;

while outline ~= -1;
    outline = fgets(outputfile);
    oline = oline + 1;
    %strncmpi(outline,'## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei',20)
    if  strncmpi(outline,'## Synthetic spectrum (Best Model) ##l_obs f_obs f_syn wei',20) == 1;
        header = oline;
    end
%    if header ~= 0 && isempty(outline) == 1;
%        nolines = oline - header;
%    end
end

%nolines
fclose('all');      
end