function bestIndex = compareIndex(f_sound,fs,f)
    f_index_1 = round(f_sound/(fs/2)*numel(f));
    f_index_2 = round(f_sound/(fs/2)*numel(f))+1;
    f_index_3 = round(f_sound/(fs/2)*numel(f))-1;
    
    a = abs(f(f_index_1) - f_sound);
    b = abs(f(f_index_2) - f_sound);
    c = abs(f(f_index_3) - f_sound);
    
    if min([a b c]) == a
        bestIndex = f_index_1;
    elseif min([a b c]) == b
        bestIndex = f_index_2;
    else
        bestIndex = f_index_3;
    end
end

