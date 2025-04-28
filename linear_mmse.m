clear;
clc;
close all;

M = 4;
data = randi([0 M-1],[2,10]);
disp('data =');
disp(data);
x = qammod(data,M,UnitAveragePower=true); % unit average power
N0 = 1;
snr = 1/N0;

disp('x');
disp(x);

Mt = 2;
Mr = 2;

H = normrnd(0, 1, [Mr,Mt]) + 1j*normrnd(0, 1, [Mr, Mt]);

disp('H');
disp(H);
fprintf('rank of H: %f\n', rank(H));

n = zeros(size(x));

y = H*x + n;
disp('y =');
disp(y);

x_hat = H'*inv((H'*H+(N0*eye(2))))*y;

rx = qamdemod(x_hat, M);
disp('rx =');
disp(rx);