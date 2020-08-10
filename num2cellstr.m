function carray = num2cellstr(numarray, formatSpec)

carray = cell(size(numarray));

for s1 = 1 : size(numarray,1)
    for s2 = 1 : size(numarray,2)
        try
            carray{s1,s2} = num2str(numarray(s1,s2), formatSpec);
        catch
            carray{s1,s2} = num2str(numarray(s1,s2));
        end
    end
end