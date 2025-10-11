
function compareStructs(a, b)
    fieldsA = fieldnames(a);
    fieldsB = fieldnames(b);
    allFields = union(fieldsA, fieldsB);

    for j = 1:length(allFields)
        f = allFields{j};
        if ~isfield(a, f)
            fprintf('  Field "%s" only in struct2\n', f);
        elseif ~isfield(b, f)
            fprintf('  Field "%s" only in struct1\n', f);
        else
            v1 = a.(f);
            v2 = b.(f);
            if ~isequal(v1, v2)
                fprintf('  Field "%s" differs:\n', f);
                disp('    struct1:'); disp(v1);
                disp('    struct2:'); disp(v2);
            end
        end
    end
end
