a = 5.43 * 10^(-10);
m_0 = 9.1 * 10^(-31);
h = 6.63 * 10^(-34);
hbar = h/(2*pi);

x = -a/2:a/100:a/2;
V = @(x) 10^(-19) * cos(2*pi/a*x) .* (x < a/4 & x > -a/4);

N = 8;

% find the fourier series coefficients of V
V_n = zeros(1,2*N + 1);

for n = -N:1:N
    V_n(n + N + 1) = real(1/a * integral(@(x) exp(-1i*2*pi*n*x/a) .* V(x), -a/2, a/2));
end

% plot the brillouin zone
k_vec = 0:pi/a:pi/a;
diag_vec = zeros(1, N + 1);
base_mat = zeros(N + 1,N + 1);
eigen_vec = zeros(N + 1, length(k_vec));

for n = -N:1:N
    series_vec = ones(1, N + 1 - abs(n)) * V_n(n + N + 1);
    base_mat = base_mat + diag(series_vec, -n);
end

for k = 1:length(k_vec)
    diag_vec = hbar ^ 2 * (k_vec(k) + (-N/2:1:N/2)*2*pi/a).^2 / (2*m_0);
    mat = base_mat + diag(diag_vec);
    eigen_vec(:,k) = eig(mat);
end

middle_gap = eigen_vec(:,1);
side_gap = eigen_vec(:,2);
middle_gap_diff = diff(middle_gap);
side_gap_diff = diff(side_gap);
gap_vec = zeros(1, N);

for i = 1:N
    if (mod(i, 2) == 1)
        gap_vec(i) = side_gap_diff(i);
    else
        gap_vec(i) = middle_gap_diff(i);
    end
end

disp(gap_vec)

estimate_gap_vec = 2 * abs(V_n(N + 2:2*N+1));

disp(estimate_gap_vec)

stem(1:N, gap_vec, 'filled')
hold on
stem(1:N, estimate_gap_vec, 'filled')
hold off

xlabel('n')
ylabel('Energy Gap (J)')
legend('Calculated Energy Gap', 'Estimated Energy Gap')

