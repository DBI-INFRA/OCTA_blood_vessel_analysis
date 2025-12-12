function res = get_struct_by_name(s, name)
    [tf, loc] = ismember(cellstr(name), {s.name});
    assert(all(tf)),'Some images were not found in previous results';
    res = s(loc);