% Error and Convergence computation with dicretization in chordwise and
% spanwise directions

close all
clear all
clc

%% Test Case 1

U_Inf_Mag =231.5;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;
Cl =[];
Cd = [];
toll = 10^(-3);
errvect_chord = [];
errvect_span = [];
errvect_span = [];
config.NCorpi = 1;
Cl2=[];
NumeroTot_vect=[];

DiscrStepChord_vect = 1:1:20;
DiscrStepSemiSpan_vect = 1:1:20;
NumeroPannelli=[];
err=1e3;

for itChord = 1:length(DiscrStepChord_vect)
    
config.RootChord = [7.88]; %[m]
config.DihedralAngle = [6]; % [°]
config.SweepAngle = [28.52]; % [°]
config.TaperRatio = [0.159]; %[-]
config.AspectRatio = [9.45]; %[-]
config.Span = [34.32]; %[m]
config.LEPosition_X = [0];  %[m]
config.LEPosition_Y = [0];  %[m]
config.LEPosition_Z = [0];  %[m]

config.RotationAngle_X = [0];
config.RotationAngle_Y = [2];
config.RotationAngle_Z = [0];

% Discretization options
config.SemiSpanwiseDiscr = [DiscrStepChord_vect(itChord)];
config.ChordwiseDiscr = [DiscrStepChord_vect(itChord)];


%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

%% Create the geometry structure

ControlPoints = cell(config.NCorpi, 1);
InducedPoints = cell(config.NCorpi, 1);
Normals = cell(config.NCorpi, 1);
InfiniteVortices = cell(config.NCorpi, 1);
Vortices = cell(config.NCorpi, 1);
internalMesh = cell(config.NCorpi, 1);
WingExtremes = cell(config.NCorpi, 1);


for iCorpo = 1:config.NCorpi

    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo);

end
    

%% Matrices initialization

NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matriceA = zeros(NPanelsTot, NPanelsTot);
matriceAA = zeros(NPanelsTot, NPanelsTot);
TermineNoto = zeros(NPanelsTot, 1);

%% Construction of the matrix
x = 0;
rowIndex = 0;

for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
   
            columnIndex = 0;
            
            ControlPointHere = ControlPoints{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            
            for jCorpo = 1:config.NCorpi

          
           
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Cycle on all of its spanwise panels
                      
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        
                        % Update column index
                        columnIndex = columnIndex + 1;
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);



                        matriceA(rowIndex, columnIndex) = dot(U, NormalHere);  
                        
                    end

                end

            end
        end
    end
end

%% Costruzione del termine noto
rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
  
            NormalHere = Normals{iCorpo}{ChordPanel_i, SpanPanel_i}.Coords;
            
            TermineNoto(rowIndex) = -dot(U_Inf, NormalHere);
            
        end
    end
end

%% Solve the linear system

Solution = linsolve(matriceA, TermineNoto);


Gamma = cell(config.NCorpi, 1);

rowIndex = 0;
for iCorpo = 1:config.NCorpi
    
    Gamma{iCorpo} = zeros( config.ChordwiseDiscr(iCorpo), config.SemiSpanwiseDiscr(iCorpo)*2 );
    
     % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iCorpo)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iCorpo)
            
            % Update row index
            rowIndex = rowIndex + 1;
            
            Gamma{iCorpo}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
        end
        
    end
    
end

%% Lift Computation
Lift2D = cell(config.NCorpi,1);


for iCorpo = 1:config.NCorpi
    L2d=[];
for Spanpanel = 1:2*config.SemiSpanwiseDiscr(iCorpo)
        l_2d= sum(Gamma{iCorpo}(:,Spanpanel)*rho*U_Inf_Mag*cosd(config.DihedralAngle(iCorpo)));
    L2d=[L2d ;l_2d];
end
Lift2D{iCorpo} = L2d;
end

%Total Lift

TotalLift= 0;
for iCorpo = 1:config.NCorpi
    L2d = Lift2D{iCorpo};
    L=sum(L2d);
    TotalLift =TotalLift + (config.SemiSpan(iCorpo)/config.SemiSpanwiseDiscr(iCorpo))*L;
end
cl = TotalLift/(0.5*rho*(U_Inf_Mag^2)*((config.Span(1)^2)/config.AspectRatio(1)));
Cl = [Cl cl];

if itChord > 1
    err = abs(Cl(itChord)-Cl(itChord-1))/abs(Cl(itChord));
    errvect_chord = [errvect_chord;err*100];
end

% Number of panels in chordwise direction
NumeroPannelli(itChord,1)=config.ChordwiseDiscr;
% Number of panels in spanwise direction
NumeroPannelli(itChord,2)=config.SemiSpanwiseDiscr*2;

NumTotPannelli=(config.ChordwiseDiscr*config.SemiSpanwiseDiscr)*2;
NumeroTot_vect=[NumeroTot_vect; NumTotPannelli];

if err <= toll
    NumTotPannelli=(config.ChordwiseDiscr*config.SemiSpanwiseDiscr)*2;
    NumeroTot_vect=[NumeroTot_vect; NumTotPannelli];
    fprintf('La convergenza per il Boeing è raggiunta con un numero di pannelli %i \n',NumTotPannelli)
    fprintf('Numero di pannelli in corda : %i\n\n',config.ChordwiseDiscr)
    fprintf('Numero di pannelli in apertura : %i\n\n',2*config.SemiSpanwiseDiscr)
    break;
end


end

figure
grid on
plot(DiscrStepChord_vect(2:itChord),errvect_chord,'*-')
title('Errore percentuale al variare del numero dei pannelli in corda Boeing 737 900ER')
xlabel('Numero pannelli in corda')
ylabel('Errore percentuale')

figure
grid on
plot(NumeroTot_vect(2:itChord),errvect_chord,'*-')
title('Errore percentuale al variare del numero dei pannelli totali Boeing 737 900ER')
xlabel('Numero pannelli totali')
ylabel('Errore percentuale')
