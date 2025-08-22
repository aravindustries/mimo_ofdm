clear;
clc;
close all;

% Same setup used in previous 2 scripts
N_FFT = 64;
N_CP = 16;
N_DATA = 48;
N_PILOT = 4;
N_NULL = 12;
N_TX = 2;
N_RX = 2;
M = 4;
pilot_indices = [7, 21, 43, 57];
data_indices = [2:6, 8:20, 22:26, 38:42, 44:56, 58:63, 1];
null_indices = [27:37, 64];
total_symbols = 100;
total_bits_per_tx = N_DATA * log2(M) * total_symbols;
snr_db = linspace(0, 50, 10);
snr = 10.^(snr_db/10);
N0 = 1./snr;

% Chan 1: Mild multipath
h11_1 = [1 0.2j 0.05-0.01j];               
h12_1 = [0.8 0.15j 0.04-0.01j];            
h21_1 = [0.7 0.1j 0.03-0.02j];             
h22_1 = [0.9 0.25j 0.06-0.01j]; 

% Chan 2: medium multipath
h11_2 = [1 0.5+0.5j 0.2 1 0.4];           
h12_2 = [0.8 0.4+0.3j 0.1 0.7 0.2];         
h21_2 = [0.7 0.3+0.4j 0.15 0.8 0.3];        
h22_2 = [0.9 0.6+0.2j 0.25 0.9 0.35];       

% Chan 3: severe!
h11_3 = [1 0.9+0.3j 0.7-0.2j 0.5+0.1j];    
h12_3 = [0.8 0.8+0.4j 0.6-0.1j 0.4+0.2j];  
h21_3 = [0.7 0.7+0.5j 0.5-0.3j 0.3+0.3j];  
h22_3 = [0.9 0.8+0.2j 0.6-0.4j 0.5+0.1j];  


% Combine channels into cells for easy manipulation
channels = {
    {h11_1, h12_1, h21_1, h22_1},
    {h11_2, h12_2, h21_2, h22_2},
    {h11_3, h12_3, h21_3, h22_3}
};


% THis step is key. We need to normalize the power of all of the channel
% gain matrixes, otherwise the channels with more taps will have a higher
% received power
for i = 1:length(channels)
    total_power = 0;
    for j = 1:4
        total_power = total_power + norm(channels{i}{j})^2;
    end
    scale = sqrt(total_power);
    for j = 1:4
        channels{i}{j} = channels{i}{j} / scale;
    end
end


channel_names = {'Mild Multipath', 'Moderate Multipath', 'Severe Multipath'};

ber_mmse = zeros(length(channels), length(snr_db));

% gen random data for both transmit antennas
data_tx1 = randi([0 M-1], 1, N_DATA * total_symbols);
data_tx2 = randi([0 M-1], 1, N_DATA * total_symbols);
data = [data_tx1; data_tx2];

% iterate thru each channel and SNR point
for ch_idx = 1:length(channels)
    h = channels{ch_idx};  % 4 channel impulse responses (each path in 2x2 MIMO)
    
    fprintf('Channel %d (%s):\n', ch_idx, channel_names{ch_idx});
    
    for snr_idx = 1:length(snr_db)
        rx_mmse = mimo_ofdm_mmse(data, M, N_FFT, N_CP, h, N0(snr_idx), pilot_indices, data_indices, N_TX, N_RX);
        errors = sum(sum(data ~= rx_mmse));
        total_bits = numel(data);
        ber_mmse(ch_idx, snr_idx) = errors / total_bits;
        
        fprintf('  SNR = %d dB: BER_MMSE = %e\n', snr_db(snr_idx), ber_mmse(ch_idx, snr_idx));
    end
end

% BER vs SNR for all channel models
figure;
for ch_idx = 1:length(channels)
    semilogy(snr_db, ber_mmse(ch_idx, :), 'o-', 'LineWidth', 2, 'DisplayName', channel_names{ch_idx});
    hold on;
end
hold off;
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title(sprintf('BER Performance of 2×2 MIMO-OFDM with MMSE Equalization M = %d', M));
legend show;


% THROUGHPUT CALCULATION
fs = 20e6; % we assume 20 MHz sampling rate

T_symbol = (N_FFT + N_CP) / fs;
bits_per_symbol = N_DATA * log2(M) * N_TX;
throughput_raw = bits_per_symbol / T_symbol;
throughput_eff = throughput_raw * (1 - ber_mmse);
disp('Effective Throughput Matrix (in Mbps):');
disp(throughput_eff / 1e6);

figure;
for ch_idx = 1:length(channels)
    plot(snr_db, throughput_eff(ch_idx, :) / 1e6, 'o-', 'LineWidth', 2, 'DisplayName', channel_names{ch_idx});
    hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('Effective Throughput (Mbps)');
title(sprintf('Effective Throughput of 2×2 MIMO-OFDM with MMSE Equalization M = %d', M));
legend show;

% MMSE Equalizer
function rx_data = mimo_ofdm_mmse(data, M, N_FFT, N_CP, h_set, N0, pilot_indices, data_indices, N_TX, N_RX)
    % reshape data for both OFDM symbols and tx antennas
    data_per_symbol = length(data_indices);
    num_symbols = size(data, 2) / data_per_symbol;
    
    % Reshape data for each tx antenna
    data_tx = cell(N_TX, 1);
    for tx = 1:N_TX
        data_tx{tx} = reshape(data(tx, :), data_per_symbol, num_symbols);
    end
    rx_data_all = cell(N_TX, 1);
    for tx = 1:N_TX
        rx_data_all{tx} = [];
    end
    
    % process all OFDM symbols
    for sym_idx = 1:num_symbols
        % QAM modulation + OFDM symbol for each tx antenna
        X_freq = zeros(N_FFT, N_TX);
        
        for tx = 1:N_TX
            % QAM modulate the data for this tx antenna
            X_data = qammod(data_tx{tx}(:, sym_idx), M, 'UnitAveragePower', true);
            
            % Create OFDM symbol
            X = zeros(N_FFT, 1);
            X(data_indices) = X_data;
            
            % Add pilots
            pilot_values = [1; 1; 1; -1];
            X(pilot_indices) = pilot_values;
            
            X_freq(:, tx) = X;
        end
        
        % IFFT to convert to time domain for each TX antenna
        x_time = zeros(N_FFT, N_TX);
        for tx = 1:N_TX
            x_time(:, tx) = ifft(X_freq(:, tx), N_FFT);
        end
        
        % CP for each TX antenna
        x_cp = [x_time(end-N_CP+1:end, :); x_time];
        
        rx_signal = zeros(length(x_cp), N_RX);
        
        % apply contributions from all tx at each rx
        for rx = 1:N_RX
            for tx = 1:N_TX
                h = h_set{(rx-1)*N_TX + tx}; 
                rx_with_h = conv(h, x_cp(:, tx));
                rx_signal(:, rx) = rx_signal(:, rx) + rx_with_h(1:length(x_cp));
            end
            
            % Add AWGN after applying channel
            noise = sqrt(N0/2) * (randn(size(rx_signal(:, rx))) + 1j*randn(size(rx_signal(:, rx))));
            rx_signal(:, rx) = rx_signal(:, rx) + noise;
        end
        
        % remove CP
        rx_no_cp = rx_signal(N_CP+1:end, :);
        
        % Convert to back to freq domain
        Y = zeros(N_FFT, N_RX);
        for rx = 1:N_RX
            Y(:, rx) = fft(rx_no_cp(:, rx), N_FFT);
        end
        
        % MMSE equalization for each subcarrier
        X_eq = zeros(N_FFT, N_TX);
        for k = 1:N_FFT
            % Skip nulls
            if ~ismember(k, data_indices) && ~ismember(k, pilot_indices)
                continue;
            end
            
            % channel for this subcarrier
            H_k = zeros(N_RX, N_TX);        
            for rx = 1:N_RX
                for tx = 1:N_TX
                    h = h_set{(rx-1)*N_TX + tx};
                    H_freq = fft([h zeros(1, N_FFT-length(h))], N_FFT);
                    H_k(rx, tx) = H_freq(k);
                end
            end
            
            % received signals for this subcarrier
            Y_k = Y(k, :).';
            
            % MMSE equalization for this subcarrier
            H_mmse = (H_k' * H_k + N0 * eye(N_TX)) \ H_k';
            X_eq_k = H_mmse * Y_k;
            
            X_eq(k, :) = X_eq_k.';
        end
        
        % demodulate data for each tx stream
        for tx = 1:N_TX
            Y_data = X_eq(data_indices, tx);
            rx_data_sym = qamdemod(Y_data, M, 'UnitAveragePower', true);
            % join the demodulated data
            rx_data_all{tx} = [rx_data_all{tx}; rx_data_sym];
        end
    end
    rx_data = zeros(size(data));
    for tx = 1:N_TX
        rx_data(tx, :) = rx_data_all{tx}.';
    end
end

