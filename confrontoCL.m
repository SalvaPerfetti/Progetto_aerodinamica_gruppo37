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
CLXfoil = [-0.2416, -0.2174, -0.1933,-0.1691 ,-0.1449, -0.1208  ,-0.0966   ,-0.0724  , -0.0483 ,-0.0241  , 0.0001, 0.0243   ,0.0484   , 0.0726   ,0.0968   ,0.1210   , 0.1451   , 0.1693   , 0.1935,  0.2176 ,0.2418 ];
% Grafico
figure;
hold on;

% Plot della formula teorica
plot(alpha_deg, Cl_theoretical, 'k-o', 'LineWidth', 2, 'DisplayName', 'Formula teorica');

% Plot dei dati di Hess-Smith
plot(alpha_HS, CLHS, 'r--', 'LineWidth', 1, 'DisplayName', 'Hess-Smith');

% Plot dei dati di Xfoil
plot(alpha_HS, CLXfoil, 'b-', 'LineWidth', 1, 'DisplayName', 'Xfoil');

% Personalizzazione del grafico
grid on;
xlabel('\alpha [deg]');
ylabel('C_l');
title('Confronto dei coefficienti di portanza (C_l)');
legend('Location', 'best');
hold off;


