% DESCRIPTION: 
% This function simulates the dynamic thermal behaviour of one heating zone in a building.
% The air, the walls and the heating system (including the production & control) are modelled
% as three lumped capacitances. The derivative is approximated by a difference quotient.
%
% INPUTS: 
%   1x11 'climate': weather inputs
%           T_ex(°C)        Exterior temperature
%           T_exgem(°C)     24h average exterior temperature
%           I_north(W/m²)   Radiation on vertical surface north
%           I_east(W/m²)    Radiation on vertical surface east
%           I_south(W/m²)   Radiation on vertical surface south
%           I_west(W/m²)    Radiation on vertical surface west
%           I_45north(W/m²) Radiation on tilted surface north
%           I_45east(W/m²)  Radiation on tilted surface east
%           I_45south(W/m²) Radiation on tilted surface south
%           I_45west(W/m²)  Radiation on tilted surface west
%           I_horiz(W/m²)   Radiation on horizontal surface
%   1x2 'occupant': internal set points and heat gains
%           T_setcomf(°C)   Setpoint for comfort temperature
%           Q_g(kJ)         Internal heat gains during time step
%   1x2 'system': heating system variables
%           Q_g(kJ)         heat gains from other systems in zone (eg SHW)
%           T_setctr(°C)    Setpoint for control sytem
%   1x3 'zoneold': previous temperature values for solving differential equation
%           T_air(°C)       air temperature at previous timestep
%           T_opaq(°C)      opaque envelope temperature at previous timestep
%           T_em(°C)        emitter temperature at previous timestep
%           T_walli(°C)     interior wall temperature at previous timestep
% PARAMETERS:
%   1x24 'building': building parameters
%           V(m³)           protected volume zone (outside measurements)
%           V_air(m³)        air volume zone (inside measurements)
%           n_50(1/h)        blowerdoor test air tightness (m³/hm²)
%           A_floor(m²)     total useful floor surface (including storeys)
%           A_floorex(m²)   total floor surface on ground or above unheated
%                           cellar (outside measurements)
%           A_wallex(m²)    total exterior wall surface (without windows)
%           A_wallint(m²)   total interior wall surface
%           A_roof(m²)      total roof surface
%           U_floor(W/m²K)  Equivalent U-value of floor slab (including
%                           thermal resistance of ground)
%           U_wall(W/m²K)   U-value of exterior walls
%           U_roof(W/m²K)   U-value of roof
%           U_window(W/m²K) U-value of window (including frame and glass)
%           Aw_north(m²)    total window surface North
%           Aw_east(m²)     total window surface East
%           Aw_south(m²)    total window surface South
%           Aw_west(m²)     total window surface West
%           Aw_45north(m²)  total tilted window surface North
%           Aw_45east(m²)   total tilted window surface East
%           Aw_45south(m²)  total tilted window surface South
%           Aw_45west(m²)   total tilted window surface West
%           Aw_horiz(m²)      total horizontal window surface
%           glassperc(-)    glass-to-window ratio
%           g-value(-)      g-value windowglass
%           cap(0/1/2)      0=light, 1=medium, 2=heavy
%   1x4 'ventilation': ventilation parameters
%           n_vent(1/h)     ventilation rate  
%           HX_eff(-)       efficiency of heat exchanger
%           n_night(1/h)    extra ventilation rate ventilative cooling (eg windows)           
%           T_night(°C)     indoor temperature at which vent. cooling is
%                           needed
%   1x9 'heating': heating system parameters
%           A_em(m²)        emitter surface (0 if not known)
%           type_em(0/1/2/3)  0=convector,1=radiator,2=floorheating,3=TABS
%           T_winregime(°C) regime temperature for water going in emitter
%           T_woutregime(°C)regime temperature for water coming out emitter
%           T_roomregime(°C)regime temperature of emitter room
%           F_reheat(W/m²)  reheat factor (16 W/m² default)
%           P_emitter(kW)   Power heat emitter (0 if not known)
%           P_prod(kW)      power heat production (0 if not known)
%           F_mod(-)        modulation ratio
%   1x3 'time'
%           step(s)         simulation timestep for building simulation
%           init(0/1)       initialisation of simulation : time = 0 => 0
%           summer(0/1)     if summer time => 1
%
% OUTPUTS: 
%   1x10 'output': outputs of model
%           T_air(°C)       air temperature
%           T_wall(°C)      wall temperature 
%           T_em(°C)        emitter temperature
%           T_op(°C)        operative temperature
%           Q_heat(kJ)      Energy transferred to emitter during timestep
%           Q_cond(kJ)      Energy lost by conduction during ts
%           Q_vent(kJ)      Energy lost by ventilation during ts
%           Q_sol(kJ)       Solar heat gains during ts
%           Q_intern(kJ)    Internal heat gains during ts
%           mcDT(kJ)        energy absorbed during ts
% 
% jeroen.van.der.veken@bbri.be; 02/06/2016


function [output]=onezone(climate,occupant,system,zoneold,building,ventilation,heating,timestep)
%% check and correct function inputs

%climate
if length(fieldnames(climate))~=11 
    disp(['Error: ',num2str(length(fieldnames(climate))),' climate args given, but 11 needed.'])
end

%occupant
if length(fieldnames(occupant))~=2 
    disp(['Error: ',num2str(length(fieldnames(occupant))),' occupant args given, but 2 needed.'])
end

%system
if length(fieldnames(system))~=2 
    disp(['Error: ',num2str(length(fieldnames(system))),' system args given, but 2 needed.'])
end

if isnumeric(system.Qsysgains)==0 || length(system.Qsysgains)==0
    system.Qg=0
    disp(['Error: System gains are set to 0']);
end
    
%zoneold
if time.init ~= 0
    if length(fieldnames(zoneold))~=3
    disp(['Error: ',num2str(length(fieldnames(zoneold))),' zoneold args given, but 3 needed.'])
    end
end

%building
if length(fieldnames(building))~=24 
    disp(['Error: ',num2str(length(fieldnames(building)),' building args given, but 24 needed.'])
end

if building.A_floor < 10 || building.A_floor > 9999
    disp(['Error: Useful floor area is not between 10m² and 9999m²'])
end

if building.V < 10 || building.V > 29999
    building.V = building.A_floor*3;
    disp(['Error: Building volume is corrected to 3 times the floor surface'])
end

if building.V_air < 0.75*building.V || building.V_air > building.V
    building.V_air = building.V*0.9;
    disp(['Error: Building air volume is corrected to 0.9 times the building volume'])
end

if building.n_50 < 0.02 || building.n_50 > 20
    building.n_50 = 7;
    disp(['Error: Building air tightness is corrected to n50=7m³/hm³'])
end

if not (any (building.cap == 0:2))
    building.cap = 1;
    disp(['Error: Building capacity type is corrected to average'])
end

%ventilation
if (length(fieldnames(ventilation))~=4 
    disp(['Error: ',num2str(length(fieldnames(ventilation))),' ventilation args given, but 4 needed.'])
end

if ventilation.n_vent < 0.5 || ventilation.n_vent > 5
    ventilation.n_vent = 0.5;
    disp(['Error: ventilation rate is corrected to minimal value of 0.5'])
end

if ventilation.HX_eff > 0.9
    ventilation.HX_eff = 0.9;
    disp(['Error: heat exchanger efficiency is corrected to maximal value of 0.9'])
end

if ventilation.HX_eff < 0
    ventilation.HX_eff = 0;
    disp(['Error: heat exchanger efficiency is corrected to minimal value of 0'])
end

%heating
if length(fieldnames(heating))~=9 
    disp(['Error: ',num2str(length(fieldnames(heating))),' heating args given, but 3 needed.'])
end

if not (any (heating.type_em == 0:3))
    heating.type_em = 1;
    disp(['Error: Emitter type is corrected to radiator'])
end

if length(heating.A_em)==0 || heating.A_em <= 0 || heating.A_em > building.A_floor
    if heating.type_em > 1
        heating.A_em=0.75*building.A_floor;
    else
        heating.A_em=building.A_floor/25;
    end
    disp(['Error: emitter surface was set to ',num2str(heating.A_em),' m²'])
end

if length(heating.F_reheat)==0 || heating.F_reheat < 0 || heating.F_reheat > 100
    heating.F_reheat=16;
    disp(['Error: reheat factor was set to ',num2str(heating.F_reheat),' W/m²'])
end

%timestep
if length(fieldnames(time))~=3 
    disp(['Error: ',num2str(length(fieldnames(time))),' timestep args given, but 3 needed.'])
end

if length(time.step)<1 || time.step<1
    time.step=60;
    disp(['Error: timestep too small and corrected to 60s'])
end

%% initialisation of simulation
if time.init < 1
    % lumped capacities
    zoneold.T_air=occupant.Tset;
    zoneold.T_opaq=occupant.Tset;
    zoneold.T_em=occupant.Tset;
    
    % building parameters
    A_opaq=building.A_floorex+building.A_wallex+building.A_roof; %m²
    UA_opaq=building.A_floorex*building.U_floor+building.A_wallex*building.U_wall+building.A_roof*building.U_roof; %W/K
    U_opaq=UA_opaq/A_opaq; %W/m²K
    R_opaqi=0.125; %m²K/W
    UA_opaqi=A_opaq/R_opaqi; %W/K
    UA_opaqe=A_opaq/(1/U_opaq-R_opaqi); %W/K
    
    A_walli=building.A_wallint; m²
    UA_walli=A_walli/R_opaqi; %W/K
    
    A_windows=building.Aw_north+building.Aw_east+building.Aw_west+building.Aw_south+building.Aw_45north+building.Aw_45east+building.Aw_45west+building.Aw_45south+building.Aw_horiz; %m²
    U_windows=building.U_window; %W/m²K
    UA_windows=A_windows*U_windows; %W/K
    
    A_floor=building.A_floor; %m²
    A_tot=A_opaq+A_windows; %m²
    UA_tot=UA_opaq+UA_windows; %W/K
    V=building.V; %m³
    V_air=building.V_air; %m³
    n_inf=building.n_50/25; %h-1
    
    C_wallex = (20+building.cap*10)*A_opaq; %kJ/K
    C_wallin = (10+building.cap*5)*(building.A_wallint+building.A_floor-building.A_floorex); %kJ/K
    C_air = 6 * V_air; %kJ/K
    C_wallexstep = C_wallex*1000/time.step; %W/K
    C_wallinstep = C_wallin*1000/time.step/1000; %W/K
    C_airstep = C_air*1000/time.step; %W/K
    
    
    % system parameters
    n_ventwinter=ventilation.n_vent*(1-ventilation.HX_eff); %h-1
    n_ventsummer=ventilation.n_vent; %h-1
    n_ventnight=ventilation.n_night; %h-1
    Hv_winter=0.34*(n_ventwinter+n_inf)*V_air; %W/K 
    
    % system dimensions
    T_exmin=-8; %°C
    P_ventmax=Hv_winter*(heating.T_roomregime-T_exmin); %W
    P_condmax=UA_tot*(heating.T_roomregime-T_exmin); %W
    P_reheat=heating.F_reheat*A_floor; %W
    P_norm=(P_ventmax+P_condmax+P_reheat)/1000; %kW
    
    if 2 > heating.P_prod || heating.Pprod < 2*P_norm %productievermogen nakijken en corrigeren ahv P_norm
        P_prod = heating.P_prod;
    else
        P_prod = 1.1*P_norm;
        disp(['Error: heat production power was reset to ',num2str(P_prod),'kW'])
    end
    
    if heating.P_emittermax < 0.8 * P_norm %emittervermogen nakijken en corrigeren ahv P_norm en P_prod
        P_emittermax = P_norm;
        disp(['Error: emitter power too small and set to ',num2str(P_emittermax),'kW'])
    elseif and(heating.P_emittermax > 1.5 * P_norm, heating.P_emittermax > P_prod)
        P_emittermax = P_norm
        disp(['Error: emitter power too large and set to ',num2str(P_prod), 'kW'])
    else P_emittermax= heating.P_emittermax
    end  
    
    if heating.type_em > 1 % surface heating
        n_em=1.1;
        if P_emittermax>100*system.A_em % max 100W/m²
            P_emittermax=100*system.A_em;
            disp(['Error: emitter power was reset to ',num2str(P_emittermax),'kW'])
        end
        if heating.type_em == 3 % TABS
           C_system = 100 * system.A_em;
           r_rad=0.55;
        else % floor heating
           C_system = 20 * system.A_em;
           r_rad=0.5;
        end
    else
        if heating.type_em == 0 % convector (0)
            C_system = 15 * P_emittermax;
            n_em=1.4;
            r_rad=0.01;
        else % radiator (1)
            C_system = 30 * P_emittermax; %kJ/K
            n_em=1.3;
            r_rad=0.15;
    end
    
    C_systemstep = C_system*1000/time.step; %W/K
    
    lndTregime = (heating.T_winregime-heating.T_woutregime) / ln((heating.T_winregime-heating.T_roomregime)/(heating.T_woutregime-heating.T_roomregime)); %K
    Cte_emheat = P_emittermax / (lndTregime^n_em); %Radiator or floor heating constant [W/K]
    R_emheat = heating.A_em * lndTregime / P_emittermax; %Resistance between emitter node and room [m²K/W]
    R_emwtosurf = R_emheat - 0.1; %Resistance between emitter node and emitter surface [m²K/W]
end

%% calculate and assign outputs

% heat gains
P_solar = building.glassperc*building.g-value*(climate.I_north*building.Aw_north+climate.I_east*building.Aw_east+climate.I_west*building.Aw_west+climate.I_south*building.Aw_south+climate.I_45north*building.Aw_45north+climate.I_45east*building.Aw_45east+climate.I_45west*building.Aw_45west+climate.I_45south*building.Aw_45south+climate.I_horiz*building.Aw_horiz); %W
P_intern = occupant.Q_g*1000/time.step; %W
P_gains = P_solar + P_intern;

% ventilation
if time.summer = 1 % if summer
    if and(zoneold.Tair > ventilation.T_night , zoneold.Tair > climate.T_ex) %conditions for free cooling
        n_venttot=n_ventsummer+n_ventnight;
    else
        n_venttot=n_ventsummer;
    end
else n_venttot=n_ventwinter;
end

Hv=0.34*n_venttot*V_air; %W/K
R_vent = 1 / Hv; %K/W
P_vent = Hv*(zoneold.Tair-climate.T_ex); %W

% conduction
P_opaqi = UA_opaqi*(zoneold.T_opaq-zoneold.T_air); %W
P_opaqe = UA_opaqe*(zoneold.T_opaq-climate.T_exgem); %W
P_windows = UA_windows*(zoneold.T_air-climate.T_ex); %W
P_walli = UA_walli*(zoneold.T_walli-zoneold.T_air);

P_losses = P_vent + P_opaqe + P_windows;

% heating


% air node
SumP_air = P_opaqi + P_walli + 0.5*P_solar + P_intern + P_emitter * (1 - r_rad) - P_vent - P_windows;









    
    
    
output.Tout=;
end