function X = create_squareobj_001(N)
X = zeros(N);

X = create_squareobj_000(N);
X(2:end-1,2:end-1) = 1-create_squareobj_000(N-2);
end