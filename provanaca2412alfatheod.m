clear
clc
x=[1.0000
    0.9500
    0.9000
    0.8000
    0.7000
    0.6000
    0.5000
    0.4000
    0.3000
    0.2500
    0.2000
    0.1500
    0.1000
    0.0750
    0.0500
    0.0250
    0.0125
         0
         0
    0.0125
    0.0250
    0.0500
    0.0750
    0.1000
    0.1500
    0.2000
    0.2500
    0.3000
    0.4000
    0.5000
    0.6000
    0.7000
    0.8000
    0.9000
    0.9500
    1.0000]
y=[ 0.0013
    0.0114
    0.0208
    0.0375
    0.0518
    0.0636
    0.0724
    0.0780
    0.0788
    0.0767
    0.0726
    0.0661
    0.0563
    0.0496
    0.0413
    0.0299
    0.0215
         0
         0
   -0.0165
   -0.0227
   -0.0301
   -0.0346
   -0.0375
   -0.0410
   -0.0423
   -0.0422
   -0.0412
   -0.0380
   -0.0334
   -0.0276
   -0.0214
   -0.0150
   -0.0082
   -0.0048
   -0.0013];
A=[x'; y'];
A=A';
dati1=A(1:18,:);
dati1= flipud(dati1);
dati2=A(19:end,:);
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
epss = eps - 0.5;

alfath = 0;
 
for i = 1:(length(epss)-1)
    alfath = alfath + (dylm(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath = (alfath*180)/pi