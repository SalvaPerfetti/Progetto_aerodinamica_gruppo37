% Cessna 172 Skyhawk

close all
clear all
clc

%% Test Case 1

U_Inf_Mag =72;
beta = 0;
U_Inf = [cosd(beta) sind(beta) 0] .* U_Inf_Mag;
rho = 1.225;
Cl =[];
Cd = [];
Cl_ala=[];
Cd_ala=[];
Cd_ala_vect=[];
Cl_ala_vect=[];
Cl_vect=[];
Cd_vect=[];

config.NCorpi = 2;
yrotvect =linspace(-10,10,21);
ycoda=[-5,0,5];

for j=1:length(ycoda)
for y = 1:length(yrotvect)

config.RootChord = [1.625 1.4];
config.DihedralAngle = [1.44 0]; % [°]
config.SweepAngle = [0 0]; % [°]
config.TaperRatio = [0.672 1]; 
config.AspectRatio = [7.52 4]; 
config.Span = [11 3.4];
config.LEPosition_X = [0 4.76];
config.LEPosition_Y = [0 0];
config.LEPosition_Z = [0 0];

config.RotationAngle_X = [0 0];
config.RotationAngle_Y = [yrotvect(y) ycoda(j)];
config.RotationAngle_Z = [0 0];

% Discretization options
config.SemiSpanwiseDiscr = [20 20];
config.ChordwiseDiscr = [20 20];


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

% Visualization
% PUNTO 1: Visualizzazione della circolazione
% figure;
% for iCorpo = 1:config.NCorpi
%     surf(1:config.SemiSpanwiseDiscr(iCorpo)*2, 1:config.ChordwiseDiscr(iCorpo), Gamma{iCorpo});
%     xlabel('Pannelli in apertura');
%     ylabel('Pannelli in corda');
%     zlabel('\Gamma (Circolazione)');
%     title('Distribuzione della circolazione sui pannelli');
%     colorbar;
%     view(3);
% end



%% Visualization of Circulation Distribution

% for iCorpo = 1:config.NCorpi
%     figure
%     [X, Y] = meshgrid(1:config.SemiSpanwiseDiscr(iCorpo)*2, 1:config.ChordwiseDiscr(iCorpo));
%     surf(X, Y, Gamma{iCorpo});
%     xlabel('Spanwise Panels');
%     ylabel('Chordwise Panels');
%     zlabel('\Gamma [Circulation]');
%     title('Circulation Distribution');
%     colormap jet;
%     colorbar;
% end



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
    if iCorpo==1
        % Within this cycle lift coefficient relative to the wing is
        % computed. Cl_ala vector contains Cl of the wing for every angle of attack with tail deflection angle fixed.
        TotalLift_ala=TotalLift;
        cl_ala=TotalLift/(0.5*rho*(U_Inf_Mag^2)*((config.Span(1)^2)/config.AspectRatio(1)));
        Cl_ala=[Cl_ala; cl_ala];
    end        
end

cl = TotalLift/(0.5*rho*(U_Inf_Mag^2)*((config.Span(1)^2)/config.AspectRatio(1)));
Cl = [Cl;cl];

% fprintf('Total Lift: %.3f N\n', TotalLift);
% Visualization of distributed lift (2D)
% figure
% for iCorpo = 1:config.NCorpi
%     plot(1:config.SemiSpanwiseDiscr(iCorpo)*2, Lift2D{iCorpo}, 'LineWidth', 1.5);
%     hold on;
% end
% xlabel('Spanwise Panels');
% ylabel('Lift [N]');
% title('Distributed Lift Along the Span');
% legend('Corpo 1');
% grid on;
% hold off;

%% quarter chord
QCcell = cell(config.NCorpi,1);
for iCorpo= 1:config.NCorpi
QC =[]; 
    for SpanPanel = 1:config.SemiSpanwiseDiscr(iCorpo)

        qcx = config.LEPosition_X(iCorpo);
        
        for Chordpanel = 1:config.ChordwiseDiscr(iCorpo)
            qqcx = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TEtip(1)-internalMesh{iCorpo}{Chordpanel,SpanPanel}.LEtip(1))/4;
            qcx = qcx + qqcx;
        end
        qcy = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TEtip(2)) + config.LEPosition_Y(iCorpo);
        qcz = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TEtip(3)) + config.LEPosition_Z(iCorpo);
        qc = [qcx qcy qcz];
        QC = [QC; qc];
    end

    SpanPanel = config.SemiSpanwiseDiscr(iCorpo);
    qxcentro = config.LEPosition_X(iCorpo);
    
    for  ChordPanel = 1:config.ChordwiseDiscr(iCorpo)
        qqxcentro = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TERoot(1)-internalMesh{iCorpo}{Chordpanel,SpanPanel}.LERoot(1))/4;
        qxcentro = qxcentro + qqxcentro;
    end
    
    qycentro = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TERoot(2)) + config.LEPosition_Y(iCorpo);
    qzcentro = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TERoot(3)) + config.LEPosition_Z(iCorpo);

    QQcentro = [qxcentro qycentro qzcentro];
    QC = [QC; QQcentro];

    for SpanPanel = config.SemiSpanwiseDiscr(iCorpo)+1:2*config.SemiSpanwiseDiscr(iCorpo)

    qcx = config.LEPosition_X(iCorpo);
        for Chordpanel = 1:config.ChordwiseDiscr(iCorpo)

            qqcx = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TEtip(1)-internalMesh{iCorpo}{Chordpanel,SpanPanel}.LEtip(1))/4;
            qcx = qcx + qqcx; 
        end
        qcy = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TEtip(2))+ config.LEPosition_Y(iCorpo);
        qcz = (internalMesh{iCorpo}{Chordpanel,SpanPanel}.TEtip(3))+ config.LEPosition_Z(iCorpo);
        qc = [qcx qcy qcz];
        QC = [QC; qc];
    end

for SpanPanel = 1: 2*config.SemiSpanwiseDiscr(iCorpo)+1
    if SpanPanel <= config.SemiSpanwiseDiscr(iCorpo)+1
    QC(SpanPanel,1) = QC(SpanPanel,1) + (config.SemiSpan(iCorpo)-(SpanPanel-1)*(config.SemiSpan(iCorpo)/config.SemiSpanwiseDiscr(iCorpo)))*tand(config.SweepAngle(iCorpo));
    else 
        QC(SpanPanel) = QC(SpanPanel) - (config.SemiSpan(iCorpo)-(SpanPanel-1)*(config.SemiSpan(iCorpo)/config.SemiSpanwiseDiscr(iCorpo)))*tand(config.SweepAngle(iCorpo));
    end
end
PM = [];
for conta=1:2*config.SemiSpanwiseDiscr(iCorpo)
    PM(conta,:)=(QC(conta,:)-QC(conta+1,:))./2+QC(conta+1,:);
end

QCcell{iCorpo}= PM;

end

Vicell=cell(config.NCorpi,1);
for iCorpo = 1:config.NCorpi
    PM = QCcell{iCorpo,1};
    Vind=[];
    
    for Spanpanel = 1:length(PM) 
        v= 0;
        ControlPointHere = PM(Spanpanel,:); 
        if iCorpo == 1
            config.NCorpi2 = 1;
        elseif iCorpo > 1
            config.NCorpi2 = config.NCorpi;
        end
        
        for jCorpo = 1:config.NCorpi2  
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo);
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        v = v + U;
                    end 
                end
              
        end 
        Vind= [Vind; v];
    end
    Vicell{iCorpo} = Vind;
end

alfaicell = cell(config.NCorpi,1);

for iCorpo = 1:config.NCorpi
    PM = QCcell{iCorpo};
    alfai = [];
    Vind = Vicell{iCorpo};
    for Spanpanel= 1:length(PM)
        Normallocal = Normals{iCorpo}{1,Spanpanel}.Coords;
        ai = atand(dot(Vind(Spanpanel,:),Normallocal)/U_Inf_Mag) + config.RotationAngle_Y(iCorpo);
        alfai = [alfai; ai];
    end
    alfaicell{iCorpo} = alfai;
end

% Computation of the induced drag

D2dcell = cell(config.NCorpi,1);
D = [];
for iCorpo = 1:config.NCorpi
    d2d = [];
    alfai = alfaicell{iCorpo};
    L2d = Lift2D{iCorpo};
    d2d(:,1) = L2d(1,:)*sind(alfai(:,1)); 
    D2dcell{iCorpo} = d2d;
    D(iCorpo) = sum(d2d);
    D_tot = (config.SemiSpan(iCorpo)/(config.SemiSpanwiseDiscr(iCorpo)))*D;
    % Computation for only the wing contribution
    if iCorpo==1
        DRAG_ala = sum(D_tot);
        cd_ala = DRAG_ala/(0.5*rho*(U_Inf_Mag^2)*((config.Span(1)^2)/config.AspectRatio(1))); 
        Cd_ala=[Cd_ala; cd_ala];
    end
end

DRAG = sum(D_tot);
cd = DRAG/(0.5*rho*(U_Inf_Mag^2)*((config.Span(1)^2)/config.AspectRatio(1)));
Cd = [Cd;cd];
end

Cd_ala_vect=[Cd_ala_vect Cd_ala];
Cl_ala_vect=[Cl_ala_vect Cl_ala];

Cl_vect=[Cl_vect Cl];
Cd_vect=[Cd_vect Cd];

Cl_ala=[];
Cd_ala=[];
Cl=[];
Cd=[];

Cd=[];
Cl=[];

end

% Plot only for the wing

    figure
for conta=1:length(ycoda)
    hold on
    grid on
    plot(yrotvect,Cl_ala_vect(:,conta))
    xlabel('\alpha (°)')
    ylabel('C_{L}')
    title('C_{L}(\alpha) (Solo Ala, Cessna 172 Skyhawk')
    legend('Deflessione = -5°','Deflessione = 0°','Deflessione = +5°')
end

    figure
for conta=1:length(ycoda)
    hold on
    grid on
    plot(yrotvect,Cd_ala_vect(:,conta))
    xlabel('\alpha (°)')
    ylabel('C_{D}')
    title('C_{D} (\alpha) (Solo Ala, Cessna 172 Skyhawk)')
    legend('Deflessione = -5° ','Deflessione = 0° ','Deflessione = +5°')
end

    figure
for conta=1:length(ycoda)
    hold on
    grid on
    plot(Cd_ala_vect(:,conta),Cl_ala_vect(:,conta))
    xlabel('C_{D}')
    ylabel('C_{L}')
    title('Polare : C_{D}(C_{L}) (Solo Ala, Cessna 172 Skyhawk)')
    legend('Deflessione = -5°','Deflessione = 0°','Deflessione = +5°')
end

% Plot for all surfaces (complete aircraft)

    figure
for conta=1:length(ycoda)
    hold on
    grid on
    plot(yrotvect,Cl_vect(:,conta))
    xlabel('\alpha (°)')
    ylabel('C_{L}')
    title('C_{L}(\alpha) (Velivolo completo, Cessna 172 Skyhawk)')
    legend('Deflessione = -5°','Deflessione = 0°','Deflessione = +5°')
end

    figure
for conta=1:length(ycoda)
    hold on
    grid on
    plot(yrotvect,Cd_vect(:,conta))
    xlabel('\alpha (°)')
    ylabel('C_{D}')
    title('C_{D}(\alpha) (Velivolo completo, Cessna 172 Skyhawk)')
    legend('Deflessione = -5°','Deflessione = 0°','Deflessione = +5°')
end

    figure
for conta=1:length(ycoda)
    hold on
    grid on
    plot(Cd_vect(:,conta),Cl_vect(:,conta))
    xlabel('C_{D}')
    ylabel('C_{L}')
    title('Polare C_{D}(C_{L}) (Velivolo completo, Cessna 172 Skyhawk)')
    legend('Deflessione = -5° ','Deflessione = 0° ','Deflessione = +5° ')
end




