function X = create_squareobj(N, obj_name, varargin)

switch upper(obj_name)
    case 'OBJ-000'
        X = create_squareobj_000(N);
    case 'OBJ-001'
        X = create_squareobj_001(N);
    case 'OBJ-002'
        X = create_squareobj_002(N);
    case 'OBJ-003'
        if length(varargin) ~= 1 
            error('To create "OBJ-003", need `k` argument'); 
        end
        k = varargin{1}; 
        X = create_squareobj_003(N, k); 
    case 'NOISE'
        X = create_squareobj_noise(N);
    otherwise
        error('"%s" is not an allowed object for creation!', obj_name);
end
end