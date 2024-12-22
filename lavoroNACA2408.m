% faccio dei tentativi con un NACA 2408 per capire quale possa essere una
% giusta approssimazione della linea media nel momento in cui aumento il
% numero di punti, noto che sia facendo una interpolazione lineare che una
% spline trovo un alfath che passa da 0.257° circa a 0.5°, dai grafici si
% nota questo a causa della pendenza iniziale nei due casi discostante dalla
% previsione usando la linea media data dalla teoria
clc
clear
close all
profilo = [
    1.00000  0.00084; 
    0.95033  0.00855;
    0.90054  0.01575;
    0.80078  0.02858;
    0.70081  0.03942;
    0.60068  0.04820;
    0.50039  0.05473;
    0.40000  0.05869;
    0.29900  0.05875;
    0.24852  0.05677;
    0.19809  0.05320;
    0.14778  0.04776;
    0.09768  0.03987;
    0.07273  0.03471;
    0.04794  0.02829;
    0.02337  0.01944;
    0.01128  0.01380;
    0.00000  0.00000;
    0.00000  0.00000;
    0.01372 -0.01134;
    0.02663 -0.01493;
    0.05206 -0.01891;
    0.07727 -0.02111;
    0.10232 -0.02237;
    0.15222 -0.02338;
    0.20191 -0.02320;
    0.25148 -0.02239;
    0.30100 -0.02125;
    0.40000 -0.01869;
    0.49961 -0.01585;
    0.59932 -0.01264;
    0.69919 -0.00942;
    0.79922 -0.00636;
    0.89946 -0.00353;
    0.94967 -0.00217;
    1.00000 -0.00084];
x_dorso = flip(profilo(1:18,1)); 
y_dorso = flip(profilo(1:18,2));
x_ventre = (profilo(19:end,1));
y_ventre = (profilo(19:end,2));

% Numero di punti da aggiungere tra ogni coppia di punti originali
n_punti = 10000; %POSSO VARIARLLA PER AGGIUNGERE DEI AIUTANO NEL CASO DELL'INTERPOLAZIONE CON UNA SPLINE

% Interpolazione lineare
x_interp_lineare = linspace(min(x_dorso), max(x_dorso), numel(x_dorso) * n_punti);
y_dorso_interp_lineare = interp1(x_dorso, y_dorso, x_interp_lineare, 'linear');
y_ventre_interp_lineare = interp1(x_ventre, y_ventre, x_interp_lineare, 'linear');

% Interpolazione spline cubica
x_interp_spline = linspace(min(x_dorso), max(x_dorso), numel(x_dorso) * n_punti);
y_dorso_interp_spline = interp1(x_dorso, y_dorso, x_interp_spline, 'spline');
y_ventre_interp_spline = interp1(x_ventre, y_ventre, x_interp_spline, 'spline');

figure;

% Dati originali
subplot(3,1, 1);
plot(x_dorso, y_dorso, 'ro-', 'DisplayName', 'Dorso originale');
hold on;
plot(x_ventre, y_ventre, 'bo-', 'DisplayName', 'Ventre originale');
grid on;
legend('Location', 'best');
title('Dati originali');
xlabel('x');
ylabel('y');

% Interpolazione lineare
subplot(3,1, 2);
plot(x_interp_lineare, y_dorso_interp_lineare, 'r-', 'DisplayName', 'Dorso (Lineare)');
hold on;
plot(x_interp_lineare, y_ventre_interp_lineare, 'b-', 'DisplayName', 'Ventre (Lineare)');
grid on;
legend('Location', 'best');
title('Interpolazione lineare');
xlabel('x');
ylabel('y');

% spline 
subplot(3,1,3);
plot(x_interp_spline, y_dorso_interp_spline, 'r-', 'DisplayName', 'Dorso (Spline)');
hold on;
plot(x_interp_spline, y_ventre_interp_spline, 'b-', 'DisplayName', 'Ventre (Spline)');
grid on;
legend('Location', 'best');
title('Interpolazione spline cubica');
xlabel('x');
ylabel('y');

%%%% Calcolo della linea media coi diversi metodi
%Linea media con interpolazione lineare
x_linea_media_lineare = x_interp_lineare; % Stessa griglia di interpolazione
y_linea_media_lineare = (y_dorso_interp_lineare + y_ventre_interp_lineare) / 2;

% Linea media con spline cubica
x_linea_media_spline = x_interp_spline; % Stessa griglia di interpolazione
y_linea_media_spline = (y_dorso_interp_spline + y_ventre_interp_spline) / 2;

% Metodo 4: Linea media teorica per un NACA 2412
m_2412 = 0.02; p_2412 = 0.4;
y_lm_naca = zeros(size(x_linea_media_lineare));
for i = 1:length(x_linea_media_lineare)
    if x_linea_media_lineare(i) <= p_2412
        y_lm_naca(i) = m_2412 / p_2412^2 * (2 * p_2412 * x_linea_media_lineare(i) - x_linea_media_lineare(i)^2);
    else
        y_lm_naca(i) = m_2412 / (1 - p_2412)^2 * ((1 - 2 * p_2412) + 2 * p_2412 * x_linea_media_lineare(i) - x_linea_media_lineare(i)^2);
    end
end

%% Plot delle linee medie calcolate
figure;

% Linea media dalle formule empiriche viste ad esercitazione
subplot(3,1,1);
plot(x_linea_media_lineare, y_lm_naca, '--', 'LineWidth', 0.8, 'DisplayName', 'Linea Media (Th)');
hold on;
plot(x_dorso, y_dorso, 'r-', 'DisplayName', 'Dorso');
plot(x_ventre, y_ventre, 'b-', 'DisplayName', 'Ventre');
grid on;
legend('Location', 'best');
title('Linea Media - DATI');
xlabel('x');
ylabel('y');

% Linea media con interpolazione lineare
subplot(3,1,2);
plot(x_linea_media_lineare, y_linea_media_lineare, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Linea Media (Lineare)');
hold on;
plot(x_dorso, y_dorso, 'r-', 'DisplayName', 'Dorso');
plot(x_ventre, y_ventre, 'b-', 'DisplayName', 'Ventre');
grid on;
legend('Location', 'best');
title('Linea Media - Interpolazione Lineare');
xlabel('x');
ylabel('y');
%lm teorica
subplot(3,1,2);
plot(x_linea_media_lineare, y_lm_naca, '--', 'LineWidth', 0.8, 'DisplayName', 'Linea Media (Th)');

% Linea media con spline cubica
subplot(3,1,3);
plot(x_linea_media_spline, y_linea_media_spline, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Linea Media (Spline)');
hold on;
plot(x_dorso, y_dorso, 'r-', 'DisplayName', 'Dorso');
plot(x_ventre, y_ventre, 'b-', 'DisplayName', 'Ventre');
grid on;
legend('Location', 'best');
title('Linea Media - Spline Cubica');
xlabel('x');
ylabel('y');
%lm teorica
subplot(3,1,3);
plot(x_linea_media_lineare, y_lm_naca, '--', 'LineWidth', 0.8, 'DisplayName', 'Linea Media (Th)');




%% Calcolo delle derivate
dydx_teoria = finiteDifferenceCentered(x_linea_media_lineare, y_lm_naca); %qui potrei avere la formula esatta ma non la uso per riuscire a validare il metodo
dydx_interplin = finiteDifferenceCentered(x_linea_media_lineare, y_linea_media_lineare);
dydx_spline =finiteDifferenceCentered(x_linea_media_spline, y_linea_media_spline);
%% Calcolo alphath
eps =x_linea_media_lineare;
epss = eps - 0.5;
alfath = 0;
for i = 1:(length(epss)-1)
    alfath = alfath + (dydx_teoria(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath_TEORIA = (alfath*180)/pi;

eps =x_linea_media_lineare;
epss = eps - 0.5;
alfath = 0;
for i = 1:(length(epss)-1)
    alfath = alfath + (dydx_interplin(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath_INTERPLIN = (alfath*180)/pi;

eps =x_linea_media_spline;
epss = eps - 0.5;
alfath = 0;
for i = 1:(length(epss)-1)
    alfath = alfath + (dydx_spline(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath_SPLINE = (alfath*180)/pi;
disp(['Alfa teoria: ', num2str(alfath_TEORIA)]); 
disp(['Alfa interpolazione lineare: ', num2str(alfath_INTERPLIN)]);
disp(['Alfa spline cubica: ', num2str(alfath_SPLINE)]);
%GLI ULTIMI DUE SI DISCOSTANO A CAUSA DELLA DERIVATA SUL BORDO D'ATTACCO
%PIU' ALTA