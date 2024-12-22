%% Hess Smith (2024)

clc
close all
clear 

addpath mat_functions

%% Input

U_inf = 1;  % Velocità all'infinito [m/s]
alpha = 2;   % Angolo di incidenza [°]
U_inf_x = U_inf * cos(deg2rad(alpha)); 
U_inf_y = U_inf * sin(deg2rad(alpha));

U_inf = [U_inf_x; U_inf_y]; % ridefinisco nelle due componenti il vettor U_inf 
U_inf_normal = [-U_inf(2); U_inf(1)]; % vettore normale alla velocità
U_inf_normal = U_inf_normal ./ norm(U_inf_normal); % versore normale alla vel

TestCase = 0; 

CodiceProfilo = '0012'; 
Chord = 1;
NPannelli = 101;

LE_X_Position = 0;
LE_Y_Position = 0;

%% Creazione profilo

% numero profilo:
Corpo = importXfoilProfile(strcat('NACA_', CodiceProfilo, '.dat')); % parte da TE e vi ritorna in senso orario
% Prima flippa i vettori
x = flipud(Corpo.x);
y = flipud(Corpo.y);
Corpo.x = x.*Chord;
Corpo.y = y.*Chord;

figure;
plot(x, y, 'o-')
axis equal

%% Creazione di una struttura di pannelli

[Centro, Normale, Tangente, Estremo_1, Estremo_2, alpha, lunghezza, L2G_TransfMatrix, G2L_TransfMatrix] = CreaStrutturaPannelli(Corpo);
        
%% Inizializzazione matrici e vettori

% Ora che ho i pannelli, posso inizializzare la matrice ed i vettori

NCols = sum(NPannelli) + 1;
NRows = NCols;
matriceA = zeros(NRows, NCols);
TermineNoto = zeros(NRows, 1);

%% Creazione della matrice quadrata As


for i = 1:NPannelli
    index_i = i; % riga

    Centro_qui = Centro(i, :)'; 
    Normale_qui = Normale(i, :)';

    indexStart_colonna = 0;

        for j = 1:NPannelli
            index_j = indexStart_colonna + j;  % Colonna
 
            Estremo_1_qui = Estremo_1(j, :)';
            Estremo_2_qui = Estremo_2(j, :)';

            L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
            G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

            matriceA(index_i, index_j) = dot(ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);

            matriceA(index_i, sum(NPannelli)+1) = matriceA(index_i, sum(NPannelli)+1) + dot(ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Normale_qui);


        end

end


%% Creazione delle componenti dei vettori a_v, c_s e c_v


Centro_Start = Centro(1, :)';
Tangente_Start = Tangente(1, :)';

Centro_End = Centro(end, :)';
Tangente_End = Tangente(end, :)';


b = 0;
for j = 1:NPannelli(1)

    index_j = j;

    Estremo_1_qui = Estremo_1(j, :)';
    Estremo_2_qui = Estremo_2(j, :)';
    L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
    G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));

    a = dot(ViSorgente(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);
    b = b + dot(ViVortice(Centro_Start, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_Start);

    a = a + dot(ViSorgente(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);
    b = b + dot(ViVortice(Centro_End, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui), Tangente_End);


    matriceA(sum(NPannelli) + 1, index_j) = a;

end

matriceA(sum(NPannelli) + 1, sum(NPannelli) + 1) = b;



%% Creazione del termine noto

for j = 1:NPannelli

    Normale_qui = Normale(j, :)';

    index = j;

    TermineNoto(index) = - dot(U_inf, Normale_qui);
end

Tangente_1 = Tangente(1, :)';
Tangente_end = Tangente(end, :)';
TermineNoto(sum(NPannelli) + 1) = - dot(U_inf, (Tangente_1 + Tangente_end));

%% Risoluzione sistema lineare
Soluzione = linsolve(matriceA,TermineNoto);
sigma = Soluzione(1:NPannelli); % Soluzione intensità sorgenti
gamma= Soluzione(NPannelli+1); % Soluzione intensità vortice

%% Calcolo del Cp e della velocità sui pannelli
V_tangenziale = zeros(NPannelli, 1);
Cp = zeros(NPannelli, 1);

for i = 1:NPannelli
    V_indotto = [0; 0];
    Centro_qui = Centro(i, :)';
    Tangente_qui = Tangente(i, :)';
    
    for j = 1:NPannelli
        Estremo_1_qui = Estremo_1(j, :)';
        Estremo_2_qui = Estremo_2(j, :)';
        L2G_TransfMatrix_qui = squeeze(L2G_TransfMatrix(j, :, :));
        G2L_TransfMatrix_qui = squeeze(G2L_TransfMatrix(j, :, :));
        
        V_indotto = V_indotto + ...
            sigma(j) * ViSorgente(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui) + ...
            gamma * ViVortice(Centro_qui, Estremo_1_qui, Estremo_2_qui, L2G_TransfMatrix_qui, G2L_TransfMatrix_qui);
    end
    
    V_tangenziale(i) = dot(U_inf, Tangente_qui) + dot(V_indotto, Tangente_qui);
    Cp(i) = 1 - (V_tangenziale(i) / norm(U_inf))^2;
end

%% Plot del coefficiente di pressione con distinzione dorso e ventre
% Identificazione dei pannelli sul dorso e sul ventre
dorso_idx = find(Centro(:, 2) > 0); % Pannelli con y > 0 (dorso)
ventre_idx = find(Centro(:, 2) <= 0); % Pannelli con y <= 0 (ventre)

figure;
hold on;
%%  Plot
% Plot del dorso
plot(Centro(dorso_idx, 1) / Chord, -Cp(dorso_idx), '-o', 'DisplayName', 'Dorso');

% Plot del ventre
plot(Centro(ventre_idx, 1) / Chord, -Cp(ventre_idx), '-o', 'DisplayName', 'Ventre');

% Impostazioni grafiche
set(gca, 'YDir', 'normal'); % Asse y normale per -Cp
xlabel('x/c');
ylabel('-C_p');
title('Distribuzione del coefficiente di pressione');
grid on;
legend('Location', 'best');
hold off;
%% Calcolo del Cl
Gamma = gamma * sum(lunghezza);
Cl = 2 * Gamma / (norm(U_inf) * Chord);
fprintf('Cl calcolato: %.4f\n', Cl);

%% Calcolo del coefficiente di momento Cm,LE

% Punto di riferimento: Bordo d'attacco (LE)
x_ref = 0; % Leading edge in x
y_ref = 0; % Leading edge in y

% Momento totale rispetto al Leading Edge
CM_LE = 0;

for i = 1:NPannelli
    % Posizione del centro del pannello rispetto al bordo d'attacco
    r_c = [Centro(i, 1) - x_ref; Centro(i, 2) - y_ref];
    
    % Normale del pannello
    n_c = Normale(i, :)';
    
    % Componente z del prodotto vettoriale (r_c x n_c)
    r_cross_n_z = r_c(1) * n_c(2) - r_c(2) * n_c(1);
    
    % Contributo del pannello al momento
    M_contrib = -Cp(i) * (lunghezza(i) / Chord^2) * r_cross_n_z;
    CM_LE = CM_LE + M_contrib;
end

% Il coefficiente di momento Cm rispetto al bordo d'attacco
Cm_LE = CM_LE;

fprintf('Il coefficiente di momento Cm,LE calcolato è: %.4f\n', Cm_LE);
%% Trasporto del coefficiente di momento al 25% della corda

% Trasporto del momento
Cm_025 = Cm_LE - Cl * 0.25;

fprintf('Il coefficiente di momento Cm,0.25 calcolato è: %.4f\n', Cm_025);
%% Cp, x da Xfoil
% Dati forniti
data = [
   1.00000  0.00126  0.40988;
   0.98964  0.00270  0.21563;
   0.97022  0.00536  0.14274;
   0.94691  0.00847  0.08773;
   0.92148  0.01177  0.04610;
   0.89522  0.01507  0.01177;
   0.86866  0.01831 -0.01732;
   0.84201  0.02146 -0.04307;
   0.81532  0.02452 -0.06670;
   0.78860  0.02749 -0.08846;
   0.76187  0.03036 -0.10934;
   0.73513  0.03314 -0.12935;
   0.70840  0.03582 -0.14884;
   0.68166  0.03840 -0.16794;
   0.65494  0.04088 -0.18708;
   0.62823  0.04325 -0.20607;
   0.60155  0.04551 -0.22520;
   0.57489  0.04764 -0.24467;
   0.54826  0.04965 -0.26422;
   0.52167  0.05152 -0.28439;
   0.49513  0.05324 -0.30473;
   0.46864  0.05481 -0.32572;
   0.44222  0.05620 -0.34699;
   0.41587  0.05740 -0.36903;
   0.38960  0.05840 -0.39152;
   0.36344  0.05918 -0.41454;
   0.33738  0.05971 -0.43835;
   0.31146  0.05999 -0.46256;
   0.28570  0.05997 -0.48759;
   0.26012  0.05964 -0.51313;
   0.23478  0.05896 -0.53957;
   0.20971  0.05790 -0.56657;
   0.18501  0.05642 -0.59428;
   0.16080  0.05449 -0.62288;
   0.13728  0.05207 -0.65230;
   0.11478  0.04914 -0.68208;
   0.09384  0.04575 -0.71185;
   0.07513  0.04203 -0.74023;
   0.05923  0.03817 -0.76491;
   0.04627  0.03439 -0.78374;
   0.03593  0.03078 -0.79453;
   0.02768  0.02739 -0.79421;
   0.02106  0.02417 -0.77917;
   0.01569  0.02107 -0.74313;
   0.01130  0.01806 -0.67766;
   0.00774  0.01507 -0.56873;
   0.00488  0.01207 -0.39955;
   0.00270  0.00905 -0.15572;
   0.00117  0.00599  0.15909;
   0.00028  0.00295  0.49751;
   0.00000  0.00000  0.77720;
   0.00028 -0.00295  0.95045;
   0.00117 -0.00599  0.99918;
   0.00270 -0.00905  0.93612;
   0.00488 -0.01207  0.80985;
   0.00774 -0.01507  0.66402;
   0.01130 -0.01806  0.52221;
   0.01569 -0.02107  0.39364;
   0.02106 -0.02417  0.28042;
   0.02768 -0.02739  0.18130;
   0.03593 -0.03078  0.09496;
   0.04627 -0.03439  0.02019;
   0.05923 -0.03817 -0.04361;
   0.07513 -0.04203 -0.09568;
   0.09384 -0.04575 -0.13664;
   0.11478 -0.04914 -0.16638;
   0.13728 -0.05207 -0.18726;
   0.16080 -0.05449 -0.20085;
   0.18501 -0.05642 -0.20896;
   0.20971 -0.05790 -0.21283;
   0.23478 -0.05896 -0.21350;
   0.26012 -0.05964 -0.21164;
   0.28570 -0.05997 -0.20779;
   0.31146 -0.05999 -0.20250;
   0.33738 -0.05971 -0.19594;
   0.36344 -0.05918 -0.18840;
   0.38960 -0.05840 -0.18019;
   0.41587 -0.05740 -0.17130;
   0.44222 -0.05620 -0.16222;
   0.46864 -0.05481 -0.15246;
   0.49513 -0.05324 -0.14261;
   0.52167 -0.05152 -0.13259;
   0.54826 -0.04965 -0.12226;
   0.57489 -0.04764 -0.11191;
   0.60155 -0.04551 -0.10114;
   0.62823 -0.04325 -0.09051;
   0.65494 -0.04088 -0.07940;
   0.68166 -0.03840 -0.06815;
   0.70840 -0.03582 -0.05642;
   0.73513 -0.03314 -0.04429;
   0.76187 -0.03036 -0.03136;
   0.78860 -0.02749 -0.01767;
   0.81532 -0.02452 -0.00264;
   0.84201 -0.02146  0.01386;
   0.86867 -0.01831  0.03261;
   0.89522 -0.01507  0.05430;
   0.92148 -0.01177  0.08115;
   0.94691 -0.00847  0.11506;
   0.97022 -0.00536  0.16152;
   0.98964 -0.00270  0.22635;
   1.00000 -0.00126  0.40988
];
%% Plot del Cp calcolato rispetto alla x/Corda
figure;
hold on;
% Plot del Cp calcolato
plot(Centro(:, 1) / Chord, -Cp, '-o', 'DisplayName', 'Cp calcolato');
xlabel('x/c');
ylabel('-C_p');
title('Distribuzione del coefficiente di pressione (calcolato)');
grid on;
legend('Location', 'best');
hold off;

%% Plot del Cp dai dati forniti rispetto alla x
% Dati forniti
x_over_c_data = data(:, 1); % x/c dai dati
Cp_data = data(:, 3);       % Cp dai dati

figure;
hold on;
% Plot del Cp dai dati forniti
plot(x_over_c_data, -Cp_data, '-o', 'DisplayName', 'Cp dai dati forniti');
xlabel('x/c');
ylabel('-C_p');
title('Distribuzione del coefficiente di pressione (dati forniti)');
grid on;
legend('Location', 'best');
hold off;

%% Grafico unico: confronto Cp calcolato e dai dati forniti
figure;
hold on;
% Cp calcolato
plot(Centro(:, 1) / Chord, -Cp, 'b-o', 'DisplayName', 'Cp calcolato con HS');
% Cp dai dati forniti
plot(x_over_c_data, -Cp_data, 'k--','LineWidth', 2, 'DisplayName', 'Cp dai dati forniti da Xfoil');
xlabel('x/c');
ylabel('-C_p');
title('Confronto Cp con \alpha= 2°');
grid on;
legend('Location', 'best');
hold off;
%% Grafico deltaCp dorso e ventre
Cp = flip(Cp); % attenzione, farlo solo 1 volta (o 3 o 5 o 2*n+1)
figure
plot(Centro(:, 1) / Chord, abs(Cp_data - Cp), 'b-o', 'LineWidth', 1.5, 'DisplayName', '$|\Delta C_p|$');

% Impostazioni del grafico
xlabel('$x/c$', 'Interpreter', 'latex');
ylabel('$|\Delta C_p|$', 'Interpreter', 'latex');
title('\textbf{Variazione del modulo di $C_p$}', 'Interpreter', 'latex'); % Titolo in grassetto
grid on;
legend('Location', 'best', 'Interpreter', 'latex');

