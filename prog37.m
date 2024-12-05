clear
clc
dati1=[1,0
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
0,0];
dati2=[0,0
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
dati1= flipud(dati1);
if dati1(:,1)==dati2(:,1)
    disp('ogay')
else
    error('no')
end

ylm = calcMidline(dati1, dati2);

plot(ylm(:,1),ylm(:,2), 'r-', 'LineWidth', 4); % Colore rosso
hold on
plot(dati1(:,1),dati1(:,2), 'b--', 'LineWidth', 2); % Colore blu tratteggiato
plot(dati2(:,1),dati2(:,2), 'g-.', 'LineWidth', 2); % Colore verde tratteggiato-punto

% Personalizzazione degli assi
axis equal; % Assicura che gli assi abbiano la stessa scala
xlabel('Asse X');
ylabel('Asse Y');
title('Grafico profilo e linea media');
legend({'linea media', 'ventre', 'dorso'}, 'Location', 'best');

% Mostra la griglia
grid on;

% Mantieni i grafici attivi
hold off;
%% CALCOLO PENDENZA LM E ANGOLO DI TH
dylm = pendenza_per_tratto(ylm(:,1),ylm(:,2));

eps = dati1(:,1);
epss = eps - 0.5

alfath = 0;
 
for i = 1:(length(epss)-1)
    alfath = alfath + (dylm(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath = (alfath*180)/pi