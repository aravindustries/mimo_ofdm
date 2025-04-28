clear;
clc;
close all;

M = 4;
data = randi([0 M-1],[2,1]);
disp('data =');
disp(data);
x_hat = qammod(data,M,UnitAveragePower=true);

disp('x_hat');
disp(x_hat);

Mt = 2;
Mr = 2;

H = normrnd(0, 1, [Mr,Mt]) + 1j*normrnd(0, 1, [Mr, Mt]);

disp('H');
disp(H);
fprintf('rank of H: %f\n', rank(H));

[U, S, V] = svd(H);

disp('Singular value Decompostion:');
disp(U);
disp(S);
disp(V);

H_hat = U*S*V';
disp('product:');
disp(H_hat);

x = V*x_hat;
disp('x =');
disp(x);

n = zeros(size(x));

y_hat = H*x + n;
disp('y_hat =');
disp(y_hat);

y = U'*y_hat;
disp('y =');
disp(y);

rx = qamdemod(y, M);
disp('rx =');
disp(rx);