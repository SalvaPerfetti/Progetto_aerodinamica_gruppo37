
clc
clear
close all
alpha = [-2.000, -1.800, -1.600, -1.400, -1.200, -1.000, -0.800, -0.600, ...
         -0.400, -0.200, 0.000, 0.200, 0.400, 0.600, 0.800, 1.000, 1.200, ...
         1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.800, 3.000, ...
         3.200, 3.400, 3.600, 3.800, 4.000, 4.200, 4.400, 4.600, 4.800, 5.000];

% Dati forniti
CL = [-0.2083, -0.1865, -0.1656, -0.1447, -0.1240, -0.1034, -0.0827, ...
      -0.0621, -0.0414, -0.0207, 0.0000, 0.0207, 0.0415, 0.0621, 0.0828, ...
       0.1034, 0.1241, 0.1448, 0.1657, 0.1865, 0.2083, 0.2304, 0.2530, ...
       0.2773, 0.3030, 0.3311, 0.3597, 0.3885, 0.4196, 0.4487, 0.4791, ...
       0.5084, 0.5383, 0.5679, 0.5976, 0.6273];

CD = [0.00696, 0.00682, 0.00669, 0.00657, 0.00647, 0.00638, 0.00631, ...
      0.00625, 0.00620, 0.00618, 0.00617, 0.00618, 0.00620, 0.00625, ...
      0.00631, 0.00638, 0.00647, 0.00657, 0.00669, 0.00682, 0.00696, ...
      0.00711, 0.00728, 0.00746, 0.00765, 0.00786, 0.00807, 0.00829, ...
      0.00853, 0.00876, 0.00901, 0.00927, 0.00954, 0.00981, 0.01010, 0.01039];

% Calcolo CL/CD e il suo massimo
CL_CD = CL ./ CD;
[max_CL_CD, idx_max] = max(CL_CD);

% Valori corrispondenti al massimo rapporto CL/CD
max_CL = CL(idx_max);
max_CD = CD(idx_max);

% Creazione del grafico
figure;
hold on;

% Plot della polare (CL vs CD)
plot(CD, CL, '-o', 'DisplayName', 'Polare (CL vs CD)');

% Retta che passa per l'origine con pendenza max(CL/CD)
x = linspace(0, max(CD)*1.2, 100);
y = (max_CL / max_CD) * x;
plot(x, y, '--r', 'DisplayName', 'Retta CL/CD max');

plot(max_CD, max_CL, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', sprintf('Max CL/CD: %.2f', max_CL_CD));

xlabel('$C_D$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$C_L$', 'Interpreter', 'latex', 'FontSize', 12);
title('Polare Aerodinamica e Retta CL/CD Massimo', 'FontSize', 14);
legend('Location', 'best');
grid on;
hold off;