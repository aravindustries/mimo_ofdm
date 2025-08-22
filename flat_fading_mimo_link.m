clear;
close all;
clc;

M = 4;
total_bits = 1000;
snr_db = linspace(0,45,15);
snr = 10.^(snr_db/10);
N0 = 1./snr; % for unit average power this is the noise power
num_iter = 1000;

% channels
H1 = [ 1,   0.14 - 0.69j;
      -0.92 + 0.17j, 0.52 - 1j ];

H2 = [ 0.5 + 1.2j,  1.1 + 0.6j;
       0.9 - 0.8j,  0.4 - 1.3j ];

H3 = [ 1 + 0j,     0.99 + 0.01j;
       0.99 - 0.01j, 1 + 0j ];

channels = {H1, H2, H3};
channel_names = {'H1', 'H2', 'H3'};

fprintf('Correlation between columns:\n');
for i = 1:3
    H = channels{i};
    h1 = H(:,1);
    h2 = H(:,2);
    rho = abs(h1' * h2) / (norm(h1) * norm(h2)); % calculate correlation
    fprintf('%s: %.4f\n', channel_names{i}, rho);
end

ber_pc_all = zeros(3, length(N0));
ber_zf_all = zeros(3, length(N0));
ber_lm_all = zeros(3, length(N0));

% iterate through each channel
for ch = 1:3
    H = channels{ch};

    ber_pc = zeros(num_iter, length(N0));
    ber_zf = zeros(num_iter, length(N0));
    ber_lm = zeros(num_iter, length(N0));

    for k = 1:num_iter
        for i = 1:length(N0)
            data = randi([0 M-1],[2,total_bits/2]);
            rx_pc = precoding(data, M, H, N0(i));
            rx_zf = zero_forcing(data, M, H, N0(i));
            rx_lm = linear_mmse(data, M, H, N0(i));

            ber_pc(k,i) = sum(data ~= rx_pc, "all") / total_bits;
            ber_zf(k,i) = sum(data ~= rx_zf, "all") / total_bits;
            ber_lm(k,i) = sum(data ~= rx_lm, "all") / total_bits;
        end
    end

    ber_pc_all(ch,:) = mean(ber_pc);
    ber_zf_all(ch,:) = mean(ber_zf);
    ber_lm_all(ch,:) = mean(ber_lm);
end

% BER PLOTS
figure;
for ch = 1:3
    subplot(3,1,ch);
    semilogy(snr_db, ber_pc_all(ch,:), '-o', 'DisplayName', 'Precoding'); hold on;
    semilogy(snr_db, ber_zf_all(ch,:), '-x', 'DisplayName', 'Zero Forcing');
    semilogy(snr_db, ber_lm_all(ch,:), '-s', 'DisplayName', 'Linear MMSE'); hold off;
    title(['BER - ', channel_names{ch}]);
    xlabel('SNR (dB)');
    ylabel('BER');
    legend();
    grid on;
end

% Throughput calculation (accounting for BER)
bits_per_symbol = log2(M);
streams = 2;
rate = bits_per_symbol * streams;

% THORUGHPUT PLOTS
figure;
for ch = 1:3
    tp_pc = rate * (1 - ber_pc_all(ch,:));
    tp_zf = rate * (1 - ber_zf_all(ch,:));
    tp_lm = rate * (1 - ber_lm_all(ch,:));

    subplot(3,1,ch);
    plot(snr_db, tp_pc, '-o', 'DisplayName', 'Precoding'); hold on;
    plot(snr_db, tp_zf, '-x', 'DisplayName', 'Zero Forcing');
    plot(snr_db, tp_lm, '-s', 'DisplayName', 'Linear MMSE'); hold off;
    title(['Throughput - ', channel_names{ch}]);
    xlabel('SNR (dB)');
    ylabel('Throughput (bits/symbol)');
    legend();
    grid on;
end

% Precoding
function rx = precoding(data, M, H, N0)
    x_hat = qammod(data,M,UnitAveragePower=true);
    [U, S, V] = svd(H); %singular value decomposition
    x = V*x_hat; % from goldsmith book
    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));
    y_hat = H*x + n; %additive white noise
    y = U'*y_hat; % from goldsmith book
    detS = det(S);
    if detS == 0
        Sinv = pinv(S);
    else
        Sinv = inv(S);
    end
    y = Sinv*y; % multiply by the inverse of the diagonal matrix
    % (intuitively makes sense, but not in the book)
    % I believe this step may be amplifying the noise in the precoding
    rx = qamdemod(y, M);
end

% zero-forcing
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
    x_hat = Hinv*y; % pretty simple and straightword
    rx = qamdemod(x_hat, M);
end

% linear-mmse
function rx = linear_mmse(data, M, H, N0)
    x = qammod(data,M,UnitAveragePower=true);
    n = sqrt(N0/2) * (randn(size(x)) + 1j*randn(size(x)));
    y = H*x + n;
    H_mmse = inv(H' * H + N0 * eye(size(H,2))) * H'; % here we account for the noise variance
    x_hat = H_mmse * y;
    rx = qamdemod(x_hat, M);
end
