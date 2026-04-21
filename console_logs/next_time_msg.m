function next_time_msg(iteration)
    persistent charCount;

    if isempty(charCount)
        charCount = 0;
    end

    fprintf(repmat('\b', 1, charCount));

    remaining = mod(iteration, 4);

    msg = "";

    switch remaining
        case 0
            msg = fprintf("( ง🔥_🔥)ง");
        case 1
            msg = fprintf("(ง🔥_🔥)ง");
        case 2
            msg = fprintf("(ง🔥_🔥) ง");
        case 3
            msg = fprintf("(ง🔥_🔥)ง");
    end

    charCount = length(msg);

end
