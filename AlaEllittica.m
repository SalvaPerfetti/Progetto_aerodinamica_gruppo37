clear;
close all;
clc;

% Input Parameters
U_Inf_Mag = 231.5;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] * U_Inf_Mag;
rho = 1.225;
Cl_alfa = 2 * pi; % Lift curve slope
alfa = 2; % degrees
lambda = 7.52; % Aspect ratio
b = 11; % Wingspan
S = (b^2) / lambda; % Wing area

% Calculations
CL_alfa = Cl_alfa / (1 + (Cl_alfa / (pi * lambda)));
alfarad = deg2rad(alfa); % Convert alfa to radians
CL = CL_alfa * alfarad;
B1 = -CL / (pi * lambda);
CD = (CL^2) / (pi * lambda);
CD_alfa = CD / alfarad;

% Initialize vectors
gammavect = [];
thetavect = 0:(pi/10):pi;


% Compute gamma for each theta
% for i = 1:length(thetavect)
%     Gamma = 2 * U_Inf_Mag * b * B1 * sin(thetavect(i));
%     gammavect = [gammavect Gamma];
% end

z_vect(:)=b*0.5*cos(thetavect(:));
Gamma_ellittica(:)=-2*U_Inf_Mag*b*B1*sqrt(1-(z_vect(:)*2/b).^2);

% Plotting gammavect vs thetavect
figure;
plot(z_vect, Gamma_ellittica(:), 'LineWidth', 1.5);
grid on
axis equal
xlabel('Posizione in apertura')
ylabel('\Gamma')
title(['Circolazione in apertura alare, ala ellittica con allungamento del ' ...
    'Cessna']);