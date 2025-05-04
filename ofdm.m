clear;
clc;
close all;

M = 16;
total_bits = 18;
data = randi([0 M-1],[1, total_bits]);

% serial to parrallel converter
N = 3;
ll_data = reshape(data, N, []);
x = qammod(ll_data,M,UnitAveragePower=true);

channel_delay_taps = 3;
u = 2;
cp = x(:, (total_bits/N-u):total_bits/N);

x = [cp, x];

X = ifft(x);





