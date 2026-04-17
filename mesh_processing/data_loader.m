function output = data_loader(filename, separator, skipheader, skipfooter)
    %DATA_LOADER Load numeric rows from a delimited text file.

    raw_data = readlines(filename);

    % Convert to string array and remove empty lines
    raw_data = string(raw_data);
    raw_data = raw_data(strlength(strtrim(raw_data)) > 0);

    n_rows = length(raw_data) - skipheader - skipfooter;

    if n_rows < 0
        error('Invalid skipheader/skipfooter for file: %s', filename);
    end

    output = cell(n_rows, 1);

    for i = 1:n_rows
        line = strtrim(raw_data(skipheader + i));

        % Split by separator and remove empty tokens
        row_strdata = strsplit(line, separator);
        row_strdata = row_strdata(~cellfun('isempty', row_strdata));

        row_numdata = zeros(1, numel(row_strdata));

        for j = 1:numel(row_strdata)
            value = str2double(row_strdata{j});

            if isnan(value)
                error('Non-numeric value in file %s, line %d.', filename, skipheader + i);
            end

            row_numdata(j) = value;
        end

        output{i} = row_numdata;
    end

end
