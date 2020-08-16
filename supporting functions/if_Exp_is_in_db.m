function [flag_equal, ix_db] = if_Exp_is_in_db( aRDS_db, db_toInclude, fields_to_check )
% In the caller function, put the following lines right after you call this
% function:
% if flag_equal == 0; continue; end;

if nargin < 3 || isempty(fields_to_check)
    fields_to_check = {};
end

flag_equal = 0;
for ii = 1 : length(db_toInclude)
    names_aRDS_db      = fieldnames(aRDS_db);
    names_db_toInclude = fieldnames(db_toInclude(ii));
    
    names_fields_int   = intersect(names_aRDS_db, names_db_toInclude);
    
    if isempty(fields_to_check)
        fields_int_to_check = names_fields_int;
    else
        fields_int_to_check   = intersect(names_fields_int, fields_to_check);
    end
    
    nr_fields = length(fields_int_to_check);
    isequal_fields = zeros(nr_fields,1);
    
    for f = 1 : nr_fields
        isequal_fields(f) = isequal(aRDS_db.(fields_int_to_check{f}), db_toInclude(ii).(fields_int_to_check{f}));
    end
    if sum(isequal_fields) == nr_fields
        flag_equal = 1;
        break
    end
%     if isequal(aRDS_db, db_toInclude(ii))
%         flag_equal = 1;
%         break
%     end
end

if flag_equal == 0
    ix_db = [];
else
    ix_db = ii;
end
