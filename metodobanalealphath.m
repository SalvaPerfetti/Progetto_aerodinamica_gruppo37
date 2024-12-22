%metodo che prende in considerazione il profilo direttamente coi dati
%assegnati e, calcolata già la linea media di quei dati si procede con
%un'integrazione
clc
clear
close all
x_corda = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.09, 0.11, 0.125, 0.15, ...
           0.18, 0.21, 0.245, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48, 0.52, 0.56, ...
           0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.83, 0.86, 0.89, 0.91, 0.925, ...
           0.94, 0.955, 0.98, 1.0];
y_lm = [0, 0.0122, 0.0166, 0.0201, 0.0231, 0.0261, 0.0318, 0.0338, 0.0355, ...
        0.0364, 0.0367, 0.036, 0.0349, 0.0333, 0.0319, 0.0303, 0.0285, ...
        0.0268, 0.0253, 0.0235, 0.0215, 0.0196, 0.0178, 0.0159, 0.0141, ...
        0.0123, 0.0106, 0.0089, 0.0076, 0.0064, 0.0051, 0.0042, 0.0035, ...
        0.0028, 0.002, 0.0009, 0]; % estratta da riga 87 di prog37

% Parametri
L = 1; % Lunghezza della corda
s = x_corda - 0.5; % passaggio di cordinate da x a s 

% Calcolo della derivata di ylm 
dylm_ds = gradient(y_lm, s); 

% Creazione della funzione di interpolazione per y'lm(s)
interp_dylm = @(s_query) interp1(s, dylm_ds, s_query, 'linear', 'extrap');

% Funzione da integrare
integrand = @(s) interp_dylm(s) ./ sqrt((L/2)^2 - s.^2);

alpha_TH = (1 / pi) * integral(integrand, -L/2, L/2);
alpha_TH =alpha_TH *180/pi;
% Risultato
disp(['Il valore di alpha_TH è: ', num2str(alpha_TH), '°']);
