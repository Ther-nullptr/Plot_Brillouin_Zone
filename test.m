mat = zeros(N + 1, N + 1);

for n = -N:1:N
    series = ones(1, N + 1 - abs(n)) * n;
    mat = mat + diag(series, -n);
end

disp(mat);