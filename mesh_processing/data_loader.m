function [output] = data_loader( ...
        filename, ...
        separator, ...
        skipheader, ...
        skipfooter ...
    )
    %   Loads numeric data from a tabular data file.
    %   The data inside the file is saved in tabular format,
    %   using the same field delimiter in each line. Each line does not need to
    %   have the same number of fields. The data are saved inside a cell array.
    %
    %   Inputs
    %   ------
    %   filename : string
    %     Name of the file with all the data.
    %   separator : string
    %     Field separator string. Typical strings are ',' for comma-separated
    %     values, '\t' for tab-separated values...
    %   skipheader : int
    %     Number of lines at the beginning of the file that do not contain
    %     numerical data to be extracted.
    %   skipfooter : int
    %     Number of lines at the end of the file that do not contain numericla
    %     data to be extracted.
    %
    %   Outputs
    %   -------
    %   output : cell(double)
    %     Cell array with all the data extracted. Each component of the array
    %     corresponds to a line of the data inside the file, converted to a
    %     one-dimensional array of numbers.

    raw_data = readlines(filename);
    output = {length(raw_data) - skipheader - skipfooter};

    for i = 1:(length(raw_data) - skipheader - skipfooter)
        row_strdata = strsplit(raw_data{skipheader + i}, separator);
        row_numdata = zeros(1, length(row_strdata));

        for j = 1:(length(row_strdata))
            row_numdata(j) = str2double(row_strdata{j});
        end

        output{i} = row_numdata;
    end

end
