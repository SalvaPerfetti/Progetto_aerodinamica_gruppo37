% metodo più semplice e attendibile per naca 24xx (ne prendo uno
% sottile NACA2408
% per coerenza con la teoria, nel senso di avere un risultato vicino il più
% possibile a quello teorico)
clear
clc
close all
M=[1.00000 0.00084;
0.95033 0.00855;
0.90054 0.01575;
0.80078 0.02858;
0.70081 0.03942;
0.60068 0.04820;
0.50039 0.05473;
0.40000 0.05869;
0.29900 0.05875;
0.24852 0.05677;
0.19809 0.05320;
0.14778 0.04776;
0.09768 0.03987;
0.07273 0.03471;
0.04794 0.02829;
0.02337 0.01944;
0.01128 0.01380;
0.00000 0.00000;
0.00000 0.00000;
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
y=M(:,2);
x=M(:,1);
A=[x'; y'];
A=A';
dati1=A(1:18,:);
dati1= flipud(dati1);
dati2=A(19:end,:);
ylm = calcMidline(dati1, dati2); %DA OSSERVARE CHE IN TAL CASO LA x/c PRESA SUL DORSO E VENTRE NON HA SIMMETRIE -> QUESTO CALCOLO è SBAGLIATO, MA NEL CASO DEL NOSTRO PROFILO E DEL NACA 0012 LE x/c son sym -> nel caso nostro va bene

plot(ylm(:,1),ylm(:,2), 'r-', 'LineWidth', 4); 
hold on
plot(dati1(:,1),dati1(:,2), 'b--', 'LineWidth', 2); 
plot(dati2(:,1),dati2(:,2), 'g-.', 'LineWidth', 2); 

axis equal; 
xlabel('Asse X');
ylabel('Asse Y');
title('Grafico profilo e linea media');
legend({'linea media', 'ventre', 'dorso'}, 'Location', 'best');

grid on;

hold off;
%% CALCOLO PENDENZA LM E ANGOLO DI TH
dylm = pendenza_per_tratto(ylm(:,1),ylm(:,2));

eps = dati1(:,1);
epss = eps - 0.5;

alfath = 0;
 
for i = 1:(length(epss)-1)
    alfath = alfath + (dylm(i)/pi)* (acos(-2*epss(i+1)) - acos(-2*epss(i)));
end
alfath = (alfath*180)/pi;