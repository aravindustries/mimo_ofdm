clear;
clc;
close all;

% 802.11a OFDM Symbol

N_FFT = 64; % FFT size (number of subcarriers!)
N_CP = 16; % CP length
N_DATA = 48; % data carriers
N_PILOT = 4; % pilots
N_NULL = 12; % null subcarriers

M = 16;

subcarrier_map = zeros(N_FFT, 1);
pilot_indices = [7, 21, 43, 57];  % pilots at -21, -7, 7, 21
data_indices = [2:6, 8:20, 22:26, 38:42, 44:56, 58:63, 1];  % data
null_indices = [27:37, 64];  % nulls

total_symbols = 100;    % num OFDM symbols
total_bits = N_DATA * log2(M) * total_symbols;
snr_db = linspace(0, 50, 10);
snr = 10.^(snr_db/10);
N0 = 1./snr;  % Noise power for unit average signal power


h1 = [1 0.01-0.005j]; % easy channel

h2 = [1 0.4+0.3j 0.15 0.2+0.1j 0.05]; % moderate multipath

h3 = [1 1.1+0.6j 0.9-0.8j 1.2+0.5j 0.8-0.3j 0.6+0.7j 0.5]; % a lot of multipath

channels = {h1, h2, h3};
channel_names = {'Less Multipath', 'Medium Multipath', 'Strong Multipath'};

% BER results arrays
ber_zf = zeros(length(channels), length(snr_db));
ber_mmse = zeros(length(channels), length(snr_db));

% random data
data = randi([0 M-1], 1, N_DATA * total_symbols);

% iterate htrough each channel and SNR point
for ch_idx = 1:length(channels)
    h = channels{ch_idx};
    
    fprintf('Channel %d (%s):\n', ch_idx, channel_names{ch_idx});
    
    for snr_idx = 1:length(snr_db)
        % Process with each equalizer

        rx_zf = ofdm_zf(data, M, N_FFT, N_CP, h, N0(snr_idx), pilot_indices, data_indices);
        ber_zf(ch_idx, snr_idx) = sum(data ~= rx_zf) / length(data);
   
        rx_mmse = ofdm_mmse(data, M, N_FFT, N_CP, h, N0(snr_idx), pilot_indices, data_indices);
        ber_mmse(ch_idx, snr_idx) = sum(data ~= rx_mmse) / length(data);
        
        fprintf('  SNR = %d dB: BER_ZF = %e, BER_MMSE = %e\n', snr_db(snr_idx), ber_zf(ch_idx, snr_idx), ber_mmse(ch_idx, snr_idx));
    end
end

% ZF vs MMSE BER
figure;
for ch_idx = 1:length(channels)
    subplot(length(channels), 1, ch_idx);
    semilogy(snr_db, ber_zf(ch_idx, :), 'g');
    hold on;
    semilogy(snr_db, ber_mmse(ch_idx, :), 'r');
    grid on;
    ylabel('BER');
    title(sprintf('Channel %d: %s', ch_idx, channel_names{ch_idx}));
    legend('Zero-Forcing', 'MMSE', 'Location', 'southwest');
    
    if ch_idx == length(channels)
        xlabel('SNR (dB)');
    end
end
sgtitle(sprintf(' OFDM ZF vs MMSE Performance by Channel for M = %d', M));
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('BER', 'FontSize', 12);

% Zero-Forcing Equalizer
function rx_data = ofdm_zf(data, M, N_FFT, N_CP, h, N0, pilot_indices, data_indices)
    data_per_symbol = length(data_indices);
    num_symbols = length(data) / data_per_symbol;
    data_reshaped = reshape(data, data_per_symbol, num_symbols);
    
    rx_data_all = [];
    
    % process each OFDM symbol
    for sym_idx = 1:num_symbols
        X_data = qammod(data_reshaped(:, sym_idx), M, 'UnitAveragePower', true);
        
        X = zeros(N_FFT, 1);
        X(data_indices) = X_data;  % place data on data subcarriers
        
        % add pilots
        pilot_values = [1; 1; 1; -1];  % 802.11a pilot pattern (BPSK)
        X(pilot_indices) = pilot_values;
        
        x = ifft(X, N_FFT); % convert to time domain
        
        % Add CP
        x_cp = [x(end-N_CP+1:end); x];
        
        % apply channel
        rx_with_cp = conv(h, x_cp);
        rx_signal = rx_with_cp(1:length(x_cp));
        
        % Add AWGN
        noise = sqrt(N0/2) * (randn(size(rx_signal)) + 1j*randn(size(rx_signal)));
        rx_signal = rx_signal + noise;
        
        % Rm cyclic prefix
        rx_no_cp = rx_signal(N_CP+1:end);
        
        % Convert back to frequency domain
        Y = fft(rx_no_cp, N_FFT);
        
        % Channel estimation (perfect CSIR)
        H = fft([h zeros(1, N_FFT-length(h))]', N_FFT);
        
        % THe magic is that the convolution with the channel in the time domain simplifies to
        % a multiplication in the frequency domain, so the channel can be
        % inverted by element-wise division when it is converted back to
        % the frequency domain!

        Y_eq = Y ./ H;
        
        % subcarriers
        Y_data = Y_eq(data_indices);
        % demod
        rx_data_sym = qamdemod(Y_data, M, 'UnitAveragePower', true);
        rx_data_all = [rx_data_all; rx_data_sym]; % append the new data
    end
    rx_data = rx_data_all';
end

%h MMSE Equalizer
function rx_data = ofdm_mmse(data, M, N_FFT, N_CP, h, N0, pilot_indices, data_indices)
    data_per_symbol = length(data_indices);
    num_symbols = length(data) / data_per_symbol;
    data_reshaped = reshape(data, data_per_symbol, num_symbols);
    rx_data_all = [];
    for sym_idx = 1:num_symbols
        X_data = qammod(data_reshaped(:, sym_idx), M, 'UnitAveragePower', true);
        X = zeros(N_FFT, 1);
        X(data_indices) = X_data;
        pilot_values = [1; 1; 1; -1];
        X(pilot_indices) = pilot_values;
        x = ifft(X, N_FFT);
        x_cp = [x(end-N_CP+1:end); x];
        rx_with_cp = conv(h, x_cp);
        rx_signal = rx_with_cp(1:length(x_cp));
        noise = sqrt(N0/2) * (randn(size(rx_signal)) + 1j*randn(size(rx_signal)));
        rx_signal = rx_signal + noise;
        rx_no_cp = rx_signal(N_CP+1:end);
        Y = fft(rx_no_cp, N_FFT);
        H = fft([h zeros(1, N_FFT-length(h))]', N_FFT);
        
        % MMSE equalizer
        H_mmse = conj(H) ./ (abs(H).^2 + N0); % noise variance is accounted for
        Y_eq = H_mmse .* Y;
        
        Y_data = Y_eq(data_indices);
        rx_data_sym = qamdemod(Y_data, M, 'UnitAveragePower', true);
        rx_data_all = [rx_data_all; rx_data_sym];
    end
    rx_data = rx_data_all';
end