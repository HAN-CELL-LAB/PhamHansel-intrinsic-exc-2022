function X = create_squareobj_002(N)
X = zeros(N);

X = create_squareobj_000(N);
k = 0;
while (N-2-k*2 > 0)
    if mod(k,2) == 0
        X(2+k:end-1-k,2+k:end-1-k) = 1-create_squareobj_000(N-2-k*2);
    else
        X(2+k:end-1-k,2+k:end-1-k) = create_squareobj_000(N-2-k*2);
    end
    k = k + 1;
end
end