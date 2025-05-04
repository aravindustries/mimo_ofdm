clear;
close all;
clc;

M = 4;
total_bits = 100;
data = randi([0 M-1],[2,total_bits/2]);

snr_db = linspace(0,30,15);
snr = 10.^(snr_db/10);
N0 = 1./snr; % for unit average power this is the noise power

H1 = [ 0.7 + 0.7j,  0.3 - 0.3j;
      -0.3 + 0.3j,  0.7 - 0.7j ];

H2 = [ 0.7 + 0.7j,  0.72 + 0.68j;
       0.7 - 0.7j,  0.68 - 0.72j ];

H3 = [ 0.9 + 0.1j,  0.85 + 0.15j;
       0.85 - 0.15j, 0.9 - 0.1j ];

rx_pc = precoding(data, M, H1, N0(6));
%rx_zf = zero_forcing(data, M, H1, N0(6));
%rx_lm = linear_mmse(data, M, H1, N0(6));

ber = sum(data ~= rx_pc, "all")/total_bits;
%disp(rx_zf);
%disp(rx_lm);

num_iter = 1000;

ber_pc = zeros(num_iter, length(N0));
ber_zf = zeros(num_iter, length(N0));
ber_lm = zeros(num_iter, length(N0));

for k=1:num_iter
    for i=1:length(N0)
        n0 = N0(i);
        data = randi([0 M-1],[2,total_bits/2]);
        rx_pc = precoding(data, M, H1, N0(i));
        rx_zf = zero_forcing(data, M, H1, N0(i));
        rx_lm = linear_mmse(data, M, H1, N0(i));
        ber_pc(k,i) = sum(data ~= rx_pc, "all")/total_bits;
        ber_zf(k,i) = sum(data ~= rx_zf, "all")/total_bits;
        ber_lm(k,i) = sum(data ~= rx_lm, "all")/total_bits;
    end
end

avg_ber_pc = mean(ber_pc);
avg_ber_zf = mean(ber_zf);
avg_ber_lm = mean(ber_lm);

figure;
semilogy(snr_db, avg_ber_pc, 'DisplayName', 'precoding');
hold on;
semilogy(snr_db, avg_ber_zf, 'DisplayName', 'zero forcing');
semilogy(snr_db, avg_ber_lm, 'DisplayName', 'linear mmse');
hold off;
legend();

function rx = precoding(data, M, H, N0)

    x_hat = qammod(data,M,UnitAveragePower=true);

    [U, S, V] = svd(H);
    x = V*x_hat;
    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));
    y_hat = H*x + n;

    y = U'*y_hat;

    detS = det(S);

    if detS == 0
        Sinv = pinv(S);
    else
        Sinv = inv(S);
    end

    y = Sinv*y;

    rx = qamdemod(y, M);
end

function rx = zero_forcing(data, M, H, N0)

    x = qammod(data,M,UnitAveragePower=true);
    
    detH = det(H);
    
    if detH == 0
        Hinv = pinv(H);
    else
        Hinv = inv(H);
    end
    
    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));

    y = H*x + n;
    x_hat = Hinv*y;
    
    rx = qamdemod(x_hat, M);
end

function rx = linear_mmse(data, M, H, N0)

    x = qammod(data,M,UnitAveragePower=true); % unit average power

    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = H*x + n;

    W_mmse = inv(H' * H + N0 * eye(size(H,2))) * H';

    x_hat = W_mmse * y;
    
    rx = qamdemod(x_hat, M);
end