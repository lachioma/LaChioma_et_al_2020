function db_int = Intersect_dbs(db1, db2, fields_to_check)

if nargin < 3 || isempty(fields_to_check)
    fields_to_check = {};
end

cnt = 0;

for i = 1 : length(db1)
    db_entry1 = db1(i);
    [flag_equal, ~] = if_Exp_is_in_db( db_entry1, db2, fields_to_check );
    if flag_equal
        cnt = cnt + 1;
        db_int(cnt,1) = db_entry1; %#ok<AGROW>
    end
end

if cnt == 0
    db_int = struct([]);
    disp('No intersection of dbs');
end
