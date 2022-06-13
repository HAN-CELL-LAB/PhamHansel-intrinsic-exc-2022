function T = filter_struct_inside_table(T, varargin)
    for i = 1:2:length(varargin)
        col_name = varargin{i};
        selected = varargin{i+1};
        T.(col_name) = filter_struct(T.(col_name), selected);
    end
end