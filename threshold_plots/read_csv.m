function structout = read_csv(filename, varargin)
% function data = read_csv(filename, ...)
%
%Reads tabular CSV data, interprets a header row, and gives you the results
%as a struct with appropriate field names. Attempts also to interpret the
%data types.
%
%Some options are: 
%
%'input': One of 'filename', 'filehandle' or 'string'. The first argument 
%is interpreted as being in this format. The default is 'filename'.
%
%'delimiter': the field delimiter character. Default ','
%
%'endofline': the line delimiter character. Default is sprintf("\n").
%
%'header':    Whether there is a header row.
%
%'skipLines': how many lines to skip before beginning.

% see also: http://wp.me/pAtBu-4O

% Peter B. Meilstrup, 2012

% This program is free software. It comes without any warranty, to
% the extent permitted by applicable law. You can redistribute it
% and/or modify it under the terms of the Creative Commons Attribution 
% 3.0 Unported (CC BY 3.0). See 
% http://creativecommons.org/licenses/by/3.0/deed.en_US for more details.

p = inputParser();
p.addRequired('filename', @ischar)
p.addParamValue('delimiter', ',', @(x)ischar(x)&&numel(x)==1);
p.addParamValue('endofline', sprintf('\n'), @(x)ischar(x)&&numel(x)==1);
p.addParamValue('input', 'filename', @(x)ismember(x,{'filename', 'filehandle', 'string'}));
p.addParamValue('output', 'struct', @(x)ismember(x,{'struct', 'dataset'}));
p.addParamValue('header', true, @islogical);
p.addParamValue('skipLines', 0, @(x) isinteger(x) && x >= 0);

p.parse(filename, varargin{:});

switch(p.Results.input)
    case 'filename'
        fid = fopen(p.Results.filename);
        onCleanup(@()fclose(fid));
        position = ftell(fid);
    case 'filehandle'
        fid = p.Results.filename;
        position = ftell(fid);
    case 'string'
        fid = p.Results.filename;
        position = 0;
end


if p.Results.skipLines > 0
    scan_with_position('%*s', 'Delimiter', p.Results.endofline, p.Results.skip);
end

if p.Results.header
    header = scan_with_position('%s', 1, 'Delimiter', p.Results.endofline);
    header = header{1};
    colnames = textscan(header{1}, '%q', 'Delimiter', p.Results.delimiter);
    colnames = colnames{1};
else
    colnames = {};
end

%textscan doesn't have any way to tell about line boundaries if you
%just use one call for the file. And I don't know a priori how many fields
%are in each line. So I have to start by reading each line, then break
%them up.
lines = scan_with_position('%s', 'Delimiter', p.Results.endofline);
lines = lines{1};
rows = cellfun( @(x)textscan(x,'%q','Delimiter',p.Results.delimiter) ...
              , lines);

n_rows = numel(rows);
n_cols = max(length(colnames), max(cellfun('prodofsize', rows)));

%populate a rectangular grid.
grid = repmat({''},n_rows, n_cols);
for r = 1:numel(lines)
    grid(r,1:numel(rows{r})) = rows{r};
end

if numel(colnames) < n_cols
    colnames = [colnames{:} repmat({''}, [1 n_cols-numel(colnames)])];
end

%Fill empty column names with spreadsheet-like letters, and ensure
%uniqueness.
letters = alphabet(n_cols);
fillwhere = cellfun(@isempty, colnames);
colnames(fillwhere) = letters(fillwhere);
colnames = genvarname(colnames);

structout = cell2struct(mat2cell(grid, n_rows, ones(1, n_cols)), colnames, 2);

%finally try converting each column to the most appropriate data type
structout = structfun(@downcast, structout, 'UniformOutput', 0);

%------------------------

function column = downcast(column, type)
    %first try boolean
    [x, what] = ismember(column, {'false', 'true', 'FALSE', 'TRUE', '0', '1', 'f', 't'});
    if all(x)
        ix = logical([0 1]);
        column = ix(mod(what, 2) + 1);
        return
    end 
    
    %then try numbers
    nums = zeros(size(column));
    pass = 1;
    for i = 1:numel(column)
        if isempty(column{i}) || strcmp(column{i}, 'NA')
            nums(i) = NaN;
        else
            [num, count, errmsg] = sscanf(column{i}, '%g');
            if (count < 1) || ~isempty(errmsg)
                nums(i) = str2double(column{i});
                if isnan(nums(i))
                    % I think this still breaks on things like '4+NaNi'
                    pass = 0;
                    break;
                end
            else
                nums(i) = num;
            end
        end
    end
    if pass
        column = nums;
    end
    %otherwise leave as strings
end

function out = scan_with_position(varargin)
    %at the beginning, we need to do consecutive reads while keeping track 
    %of the position, which textscan handles differently for strings vs. 
    %filehandles. Ugh. This is why R has string connections.
    switch(p.Results.input)
        case 'filename'
            %we would fseek but don't need to.
            out = textscan(fid, varargin{:});
        case 'filehandle'
            %we would fseek but don't need to.
            out = textscan(fid, varargin{:});
        case 'string'
            [out, nread] = textscan(fid(position+1:end), varargin{:});
            position = position + nread;
    end
end

function letters = alphabet(k)
    %generate a sequence of letters like spreadsheet columns
    
    %spreadsheet counting is a bit weird, it goes A..Z, then AA..ZZ, then
    %AAA..ZZZ, etc.
    prelude = '';
    nchars = 1;
    while 26^nchars < k
        prelude = strvcat(prelude, dec2base(0:(26^nchars - 1), 26));
        k = k - (26^nchars);
        nchars = nchars + 1;
    end

    base26 = strvcat(prelude, dec2base(0:k-1, 26, nchars));
    lettermask = ismember(base26, 'A':'Z');
    numbermask = ismember(base26, '0':'9');
    letters = base26;
    letters(numbermask) = base26(numbermask) - '0' + 'A';
    letters(lettermask) = base26(lettermask) - 'A' + 'K';
    letters = cellstr(letters);
end

end