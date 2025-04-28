clear;
clc;
close all;

M = 4;
data = randi([0 M-1],[2,1]);
disp('data =');
disp(data);
x = qammod(data,M,UnitAveragePower=true);

disp('x');
disp(x);

Mt = 2;
Mr = 2;

H = normrnd(0, 1, [Mr,Mt]) + 1j*normrnd(0, 1, [Mr, Mt]);

disp('H');
disp(H);
fprintf('rank of H: %f\n', rank(H));
detH = det(H);
fprintf('det(H) = %f\n', detH);

if detH == 0
    Hinv = pinv(H);
else
    Hinv = inv(H);
end

n = zeros(size(x));

y = H*x + n;
disp('y =');
disp(y);

x_hat = Hinv*y;

rx = qamdemod(x_hat, M);
disp('rx =');
disp(rx);