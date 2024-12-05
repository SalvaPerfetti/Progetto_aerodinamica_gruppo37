% Definizione del vettore di angoli di attacco alpha
clc
clear
close all
alpha_deg = -2:0.1:2; % Vettore da -2° a 2° con passo 0.1
alpha_rad = deg2rad(alpha_deg); % Conversione in radianti

% Calcolo dei coefficienti di portanza Cl con la formula teorica
Cl_theoretical = 2 * pi * alpha_rad; % Formula teorica per un profilo sottile

% Dati Cl da Hess-Smith
alpha_HS = -2:0.2:2; % Vettore di alpha per Hess-Smith
CLHS = [-0.2464 -0.2218 -0.1971 -0.1725 -0.1479 -0.1232 -0.0986  -0.0739 ...
        -0.0493 -0.0246 0.000 0.0246 0.0493 0.0739 0.0986 0.1232 0.1479 ...
         0.1725 0.1971 0.2218 0.2464]; 

% Dati Cl da Xfoil
CLXfoil = [-0.2083, -0.1865, -0.1656, -0.1447, -0.1240, -0.1034, -0.0827, ...
           -0.0621, -0.0414, -0.0207, 0.0000, 0.0207, 0.0415, 0.0621, ...
            0.0828, 0.1034, 0.1241, 0.1448, 0.1657, 0.1865, 0.2083];

% Grafico
figure;
hold on;

% Plot della formula teorica
plot(alpha_deg, Cl_theoretical, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Formula teorica');

% Plot dei dati di Hess-Smith
plot(alpha_HS, CLHS, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Hess-Smith');

% Plot dei dati di Xfoil
plot(alpha_HS, CLXfoil, 'k-o', 'LineWidth', 1.5, 'DisplayName', 'Xfoil');

% Personalizzazione del grafico
grid on;
xlabel('\alpha [deg]');
ylabel('C_l');
title('Confronto dei coefficienti di portanza (C_l)');
legend('Location', 'best');
hold off;


