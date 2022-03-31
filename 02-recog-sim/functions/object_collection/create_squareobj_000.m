function X = create_squareobj_000(N)
X = zeros(N);

X(1:round(N/2),1:round(N/2)) = 1;
X(round(N/2)+1:end,round(N/2)+1:end) = 1;

end