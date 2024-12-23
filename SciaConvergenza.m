close all
clear all
clc

%% Test Case 1

U_Inf_Mag =72; %m/s
beta = 0; %deg
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag; %m/s
rho = 1.225; %kg/m^3
Cl =[];
Cd = [];
config.NCorpi = 1;
yrotvect = 40:5:70;
LL = [];
toll = 10^(-3);
for n = 1:length(yrotvect)

config.RootChord = [3]; %
config.DihedralAngle = [-5]; % [°]
config.SweepAngle = [15]; % [°]
config.TaperRatio = [1/3]; %[-]
config.AspectRatio = [6.2]; %[-]
config.Span = [12]; %[m]
config.LEPosition_X = [0]; %[m]
config.LEPosition_Y = [0]; %[m]
config.LEPosition_Z = [0]; %[m]

config.RotationAngle_X = [0];%[deg]
config.RotationAngle_Y = [2];%[deg]
config.RotationAngle_Z = [0];%[deg]

% Discretization options
config.SemiSpanwiseDiscr = [20];
config.ChordwiseDiscr = [20];


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
    [ControlPoints{iCorpo}, InducedPoints{iCorpo}, Normals{iCorpo}, InfiniteVortices{iCorpo}, Vortices{iCorpo}, internalMesh{iCorpo}, WingExtremes{iCorpo}] = createStructure(config, iCorpo,yrotvect(n));
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
LL = [LL TotalLift];

end

% Compute the error and find the position of the value for which the
% problem converges

errvect =[]; 

for it = 1:length(LL)-1
    err = abs(LL(it+1)-LL(it))/abs(LL(it+1));
    errvect = [errvect err];
    if err<toll
        fprintf('La scia converge per uno lungo %i volte la misura della corda',yrotvect(it+1))
        break
    end
end

