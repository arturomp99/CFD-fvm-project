function [new_results] = sample_results(w, t, old_results, ...
        sampling_interval)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    new_results = old_results;

    if (size(old_results, 1) > 0)
        last_t = old_results(end, 1);
    else
        last_t = -1E100;
    end

    if ((t - last_t) >= sampling_interval)
        new_results = [old_results; [t, w']];
    end

end
