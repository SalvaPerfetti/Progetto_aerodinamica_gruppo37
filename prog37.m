clc
clear
close all
profilo=[1,0
0.98,-0.002
0.955,-0.0043
0.94,-0.0057
0.925,-0.007
0.91,-0.0084
0.89,-0.01
0.86,-0.0125
0.83,-0.0149
0.8,-0.0172
0.76,-0.0198
0.72,-0.0222
0.68,-0.0244
0.64,-0.0263
0.6,-0.0279
0.56,-0.0293
0.52,-0.0304
0.48,-0.031
0.44,-0.0314
0.4,-0.0315
0.36,-0.031
0.32,-0.0299
0.28,-0.0281
0.245,-0.026
0.21,-0.0231
0.18,-0.02
0.15,-0.0167
0.125,-0.0142
0.11,-0.013
0.09,-0.0115
0.075,-0.011
0.05,-0.0108
0.04,-0.0108
0.03,-0.01
0.02,-0.009
0.01,-0.0077
0,0
0,0
0.01,0.0321
0.02,0.0423
0.03,0.0503
0.04,0.0571
0.05,0.0631
0.075,0.0745
0.09,0.0791
0.11,0.084
0.125,0.0869
0.15,0.0901
0.18,0.092
0.21,0.0929
0.245,0.0926
0.28,0.092
0.32,0.0905
0.36,0.088
0.4,0.085
0.44,0.082
0.48,0.078
0.52,0.0733
0.56,0.0685
0.6,0.0635
0.64,0.0582
0.68,0.0527
0.72,0.0469
0.76,0.041
0.8,0.035
0.83,0.0301
0.86,0.0253
0.89,0.0203
0.91,0.0167
0.925,0.014
0.94,0.0113
0.955,0.0084
0.98,0.0038
1,0];
profilo=flip(profilo);
x_dorso = flip(profilo(1:37,1)); 
y_dorso = flip(profilo(1:37,2));
x_ventre = (profilo(38:end,1));
y_ventre = (profilo(38:end,2));
%calcolo l.m. dai dati
dati1= [x_dorso y_dorso];
dati2= [x_ventre y_ventre];
ylm = calcMidline(dati1, dati2);
ylm=ylm(:,2);

% Numero di punti da aggiungere tra ogni coppia di punti originali
n_punti = 10000;

% Interpolazione lineare
x_interp_lineare = linspace(min(x_dorso), max(x_dorso), numel(x_dorso) * n_punti);
y_dorso_interp_lineare = interp1(x_dorso, y_dorso, x_interp_lineare, 'linear');
y_ventre_interp_lineare = interp1(x_ventre, y_ventre, x_interp_lineare, 'linear'); % praticamente essendo una interpolazione lineare darà lo stesso grafico di quello prendendo solo i dati (anche lo stess alphath quindi), lka uso esclusivamente per ampliare il numero dei punti lungo la qualcìe posso fare la diff degli integrandi per capire dove sono i discostamenti più alti (li aspetterò sul BA e su di esso mi aspetterò anche una curva più ripida a causa del peso 1/sqrt(L^2*0.25 -s^2)-> inf)

% 
x_interp_2 = linspace(min(x_dorso), max(x_dorso), numel(x_dorso) * n_punti);
y_dorso_interp_2 = interp1(x_dorso, y_dorso, x_interp_lineare, 'makima'); % simile alla spline ma preserva le possibili oscillazioni (infatti da un alfa più basso e vicino a quello dei dati)
y_ventre_interp_2 = interp1(x_ventre, y_ventre, x_interp_lineare, 'makima');

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
plot(x_interp_2, y_dorso_interp_2, 'r-', 'DisplayName', 'Dorso Makima');
hold on;
plot(x_interp_2, y_ventre_interp_2, 'b-', 'DisplayName', 'Ventre Makima');
grid on;
legend('Location', 'best');
title('Interpolazione Makima');
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
%Linea media con interpolazione 2
x_linea_media_2 = x_interp_2; % Stessa griglia di interpolazione
y_linea_media_2 = (y_dorso_interp_2 + y_ventre_interp_2) / 2;



% Linea media con spline cubica
x_linea_media_spline = x_interp_spline; % Stessa griglia di interpolazione
y_linea_media_spline = (y_dorso_interp_spline + y_ventre_interp_spline) / 2;



%% Plot delle linee medie calcolate
figure;

% Linea media trivial
subplot(3,1,1);
plot(x_dorso, ylm, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Linea Media dati');
hold on;
plot(x_dorso, y_dorso, 'r-', 'DisplayName', 'Dorso');
plot(x_ventre, y_ventre, 'b-', 'DisplayName', 'Ventre');
grid on;
legend('Location', 'best');
title('Linea Media - Dati');
xlabel('x');
ylabel('y');

% Linea media con makima
subplot(3,1,2);
plot(x_linea_media_2, y_linea_media_2, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Linea Media Makima');
hold on;
plot(x_interp_2, y_dorso_interp_2, 'r-', 'DisplayName', 'Dorso');
plot(x_interp_spline, y_ventre_interp_spline, 'b-', 'DisplayName', 'Ventre');
grid on;
legend('Location', 'best');
title('Linea Media - Makima');
xlabel('x');
ylabel('y');

% Linea media con spline cubica
subplot(3,1,3);
plot(x_linea_media_spline, y_linea_media_spline, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Linea Media (Spline)');
hold on;
plot(x_interp_spline, y_dorso_interp_spline, 'r-', 'DisplayName', 'Dorso');
plot(x_interp_spline, y_ventre_interp_spline, 'b-', 'DisplayName', 'Ventre');
grid on;
legend('Location', 'best');
title('Linea Media - Spline Cubica');
xlabel('x');
ylabel('y');




%% Calcolo delle derivate

dydx_lm_interp1 = finiteDifferenceCentered(x_linea_media_lineare , y_linea_media_lineare );
dydx_lm = finiteDifferenceCentered(x_dorso, ylm);
dydx_interp2 = finiteDifferenceCentered(x_linea_media_2, y_linea_media_2);
dydx_spline =finiteDifferenceCentered(x_linea_media_spline, y_linea_media_spline);
%% Calcolo alphath

eps = x_linea_media_lineare ;
epss = eps - 0.5;

alfath = 0;
 
for i = 1:(length(epss)-1)
    alfath = alfath + (dydx_lm_interp1 (i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath_dati = (alfath*180)/pi;

eps =x_linea_media_2;
epss = eps - 0.5;
alfath = 0;
for i = 1:(length(epss)-1)
    alfath = alfath + (dydx_interp2(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath_INTERP2 = (alfath*180)/pi;

eps =x_linea_media_spline;
epss = eps - 0.5;
alfath = 0;
for i = 1:(length(epss)-1)
    alfath = alfath + (dydx_spline(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath_SPLINE = (alfath*180)/pi;
disp(['Alfa dati: ', num2str(alfath_dati)]);
disp(['Alfa interpolazione usando makima: ', num2str(alfath_INTERP2)]);
disp(['Alfa spline cubica: ', num2str(alfath_SPLINE)]); 
%% Primo set di grafici: Differenze tra le derivate
figure;

plot(x_linea_media_lineare, abs(dydx_lm_interp1 - dydx_interp2), 'r-', 'LineWidth', 1.5, ...
    'DisplayName', '$\left| \frac{dy_{\mathrm{interp1}}}{dx} - \frac{dy_{\mathrm{makima}}}{dx} \right|$');
hold on;
plot(x_linea_media_lineare, abs(dydx_lm_interp1 - dydx_spline), 'b--', 'LineWidth', 1.5, ...
    'DisplayName', '$\left| \frac{dy_{\mathrm{interp1}}}{dx} - \frac{dy_{\mathrm{spline}}}{dx} \right|$');
grid on;
legend('Location', 'best', 'Interpreter', 'latex');
title('Differenze tra derivate', 'Interpreter', 'latex');
xlabel('$\frac{x}{c}$', 'Interpreter', 'latex');
ylabel('$\\| \\{Differenza} \\|$', 'Interpreter', 'latex');



%% Secondo set di grafici: Integrandi della formula di Theodorsen
L = 1; % Lunghezza della corda
s = x_linea_media_lineare - 0.5; % Traslazione di x rispetto al centro della corda
integrand_lm_interp1 = dydx_lm_interp1 ./ sqrt((L^2)/4 - s.^2);
integrand_interp2 = dydx_interp2 ./ sqrt((L^2)/4 - s.^2);
integrand_spline = dydx_spline ./ sqrt((L^2)/4 - s.^2);

figure;
hold on;
plot(x_linea_media_lineare, abs(integrand_lm_interp1-integrand_interp2), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Integrando (makima)');
plot(x_linea_media_lineare, abs(integrand_lm_interp1-integrand_spline), 'k-.', 'LineWidth', 1.5, 'DisplayName', 'Integrando (spline)');
grid on;
legend('Location', 'best');
title('Differenza degli Integrandi della formula di Theodorsen');
xlabel('x/c');
ylabel('Integrando');