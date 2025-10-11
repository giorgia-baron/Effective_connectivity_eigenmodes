function compareStructArrays(s1, s2)
    if length(s1) ~= length(s2)
        error('Struct arrays are of different lengths.');
    end

    n = length(s1);
    for i = 1:n
        if ~isequal(s1(i), s2(i))
            fprintf('Structs differ at row %d:\n', i);
            pause
            compareStructs(s1(i), s2(i));  % Call field-by-field diff (defined below)
        end
    end
end
