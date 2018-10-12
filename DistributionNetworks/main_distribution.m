%% Data-based distributionally robust stochastic optimal power flow (OPF)
% Overvoltage mitigation for IEEE 37-node distribution power network
%
% Latest update: 2018-Oct-12
%
% Copyright (c) 2018-2019, Yi Guo, Tyler Summers (PI)
%
% Control, Optimization and Networks Laboratory
% Department of Mechanical Engineering
% The University of Texas at Dallas, Richardson, TX, USA.
%
% Emails:   yi.guo2@utdallas.edu,
%           tyler.summers@utdallas.edu.
%
% This script provides a computationally-affordable chance-contrained AC 
% OPF framework to solve an overvoltage problem for IEEE 37-node 
% distribution network. The proposed distributionally robust stochastic OPF 
% methodology mitigate overvoltages by controlling set points in renewable 
% energy resources and energy storage devices. The set points of 
% controllable devices are repeatedly optimized over a finite planning 
% horizon within a MPC feedback scheme. The risk conservativeness of the 
% voltage magnitude constraints and the out-of-sample performance 
% robustness to sampling errors are explicitly adjustable by two scalar 
% parameters (e.g., weight factor and Wasserstein distance).
% 
% In general, this framework works for any solar power resolutions (e.g., 30 min)
% and forecasting methods, please refer to "README.md" for how to import
% necessary data: solar forecasting and electric power loads.
%
% The mathematical formulation and results discussions can be found in the
% following two references.
%
% [1]Y.Guo, K.Baker, E.Dall'Anese, Z.Hu and T.Summers, "Data-based
% distribubtionally robust stochastic optimal power flow, Part
% I: Methodologies", available on arXiv.org.
%
% [2]Y.Guo, K.Baker, E.Dall'Anese, Z.Hu and T.Summers, "Data-based
% distribubtionally robust stochastic optimal power flow, Part
% II: Case studies", available on arXiv.org.
%
% This material is based on works supported by the National Science
% Foundation (NSF) under grant CNS-1566127. (PI: Dr. Tyler H. Summers)
%
% Questions and comments are welcome.
%
% MIT License
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

clc; clear all; close all;
% warning off

disp(['Running Data-based distributionally robust stochastic optimal power flow (OPF)'])
disp(['Overvoltage problem for IEEE 37-node distribution power networks'])
disp(['Latest update: 2018-Oct-12'])

disp(['..................................................................................'])

disp(['Copyright (c) 2018-2019, Yi Guo, Tyler Summers (PI)'])
disp(['Control, Optimization and Network Laboratory'])
disp(['Department of Mechanical Engineering'])
disp(['The University of Texas at Dallas, Richardson, TX, USA.'])

%% load input data
disp(['..................................................................................'])
load('solar.mat')
disp(['solar data load sucessfully ......'])

load('loaddata.mat')
disp(['load data load sucessfully ......'])
disp(['..................................................................................'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% program and system setup %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = 1;
Timeinterval = 5; % interval 5 mins
tot_time = 1440/Timeinterval; % total timesteps within 24 hours with 5 min timeinterval

rhoMatrix = [1,5,10,15,20,30,35,40,100,200];% formulate a vector with multiple weight factors \rho
varepsilon = 0.0000; % Wasserstein distance \espilon
horizon = 3; % 12 = 1 hour; pick a value in n_time
n_time = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 21, 24, 30, 36, 42, 48];
a_time = [5 10 15 20 25 30 40 50 60 75 90 105 120 150 180 210 240];
phase = 3; % set three-phase electric wires
ndis = 0.95; % Battery discharging efficiency
nch = 0.95; % Battery charging efficiency
Nnode = 36; % number of nodes
Nlines = 43; % number of lines
Zbase = 1; % Impedence Base is 1;
Vbase =  4800; % Voltage base
Sbase = (Vbase^2)/Zbase; % power base
MagPcc = 1.02;% Known voltage at the secondary of the transformer
theta_0 = 0; % Reference angles

epsilon = 0.01; % Constraint fulfillment parameter for voltage constraints.
nu = 0.01; % Constraint fulfillment parameter for inverter power constraint

tot_cvx_cost = 0;

% Voltage limits in p.u.
Vmin = 0.95; % voltage magnitude lower bound in p.u.
Vmax = 1.05; % voltage magnitude upper bound in p.u.

%% Loads %%
% Nodes with loads
Node_P_l = [2 5 6 7 9 10 11 13 14 19 20 21 22 24 26 27 28 29 30 ...
    32 33 35 36]; % locations of  active loads
Node_Q_l = Node_P_l; % locations of reactive loads

Nloads = length(Node_P_l); % number of loads

%% renewable energy resources (RESs) setup (PVs)

% PV locations
Node_PV = [4 7 9 10 11 13 16 17 20 22 23 26 28 29 30 31 32 33 34 35 36];

% The number of PVs
Ninverters = length(Node_PV);

% capacity of PVs in kVA
maximum_S = [150 300 300 600 660 360 600 360 450 150 750 300 750 300 360 ...
    600 330 750 450 450 450];

% Nodes without any PV systems
Node_noPV = setdiff(2:Nnode,Node_PV);

%% Electricity Costs/Feed-in Tariffs %%
% Time dependent
c_i = 10*ones(1, tot_time+horizon*2); % cost of power from utility
d_i = 3*ones(1, tot_time+horizon*2); % feed-in tarrif cost to utility
e_i = 3*ones(1, tot_time+horizon*2); % reactive power generation cost for inverters
f_i = e_i*2;

%% Battery storage setup

% locations of energy storage devices
Node_batt = [9 10 28 29 32 35 36]; % shown in the paper
Node_nobatt = setdiff(1:Nnode,Node_batt);
Nbatt = length(Node_batt); % Number of batteries

maximum_S_bat = [100 100 50 250 250 120 200]; % in kWh shown in the paper

Bmax = maximum_S_bat(1:Nbatt).'*1000/Sbase; % upper bound of state of charge
Bmin = 0.*Bmax; % lower bound of state of charge
B_0 = Bmin;% Initialize the state of charge

Delta_t = 5; % Timestep for storage decisions

%% number of samples in dataset %%
NumS = 30; % for distributionally robust optimization, N.
NumMC = 92; % The number of scenarios for Monte Carlo simulation.
%% normalize RESs output in p.u.
PV_Real24 = PV_Real24./Sbase; % PV_Real24(Timestep,1)
PV_Errors_MPC_S = PV_Errors_MPC_S./Sbase; % forecast errors PV_Errors_MPC_S(Timestep,Horizon,Scenarios)
PV_MC = PV_MC./Sbase; % solar data prepared for Monte Carlo Simulation PV_MC(Timestep,Scenarios)

PV_Real24 = [PV_Real24 PV_Real24(1:horizon)];

% Distribute solar data over Node_PV nodes
totLoad = 0;
totPV = 0;

% prepare and check solar data for monte carlo simulation
for i = 1 : Ninverters
    for j = 1 : NumMC
        PV_Aggregate_S(j,i,:) = PV_MC(:,j);
    end
end

for i = 1 : Ninverters
    for j = 1 : NumMC
        for tt = 1 : tot_time
            if (PV_Aggregate_S(j,i,tt) > maximum_S(i)./Sbase*1000);
                PV_Aggregate_S(j,i,tt) = maximum_S(i)./Sbase*1000;
            end
        end
    end
end

% prepare and check solar data at distributionally robust optimization
for i = 1 : Ninverters
    PV_Aggregate(i,:) = PV_Real24;
end

for p = 1 : size(PV_Aggregate,1)
    for k = 1 : size(PV_Aggregate,2)
        if(PV_Aggregate(p,k) > maximum_S(p)./Sbase*1000)
            PV_Aggregate(p,k) = maximum_S(p)./Sbase*1000;
        end
    end
end

PV_MC_S = PV_Aggregate_S;
PV_SDP = PV_Aggregate; % solar data for solving voltage SDP relaxation

PV_Aggregate1 = zeros(Nnode,length(PV_Real24));
PV_Aggregate1(Node_PV,:) = PV_Aggregate;
PV_Aggregate = PV_Aggregate1;

totLoad = sum(P_l,1); % marked aggregated load over feeder
totPV = sum(PV_Aggregate,1); % marked aggregate PV over feeder

% plot aggregated load and solar
figure;
plot(totLoad(1:tot_time)*Sbase/1000,'LineWidth',2)
hold on
plot(totPV(1:tot_time)*Sbase/1000,'LineWidth',2)
legend('Total Loads','Total PV')
xlabel('Time')
ylabel('kW')

% plot all solar scenarios for monte carlo simulation
figure;
for ss = 1:NumMC
    plot(sum(squeeze(PV_MC_S(ss,:,:))*Sbase/1000,1))
    hold on
    plot(totPV(1:tot_time)*Sbase/1000,'LineWidth',3)
end
xlabel('Time')
ylabel('kW')
title('All scenarios for Monte Carlo Simulation')

%% Admittance Matrix %%%%%%%%%
[Y_net] = form_admittance(phase,Zbase,Sbase);

% Create the matrices for trace(Y*V)
Yp = zeros(Nnode,Nnode*Nnode);
Yq = zeros(Nnode,Nnode*Nnode);
M = zeros(Nnode,Nnode*Nnode);
for nn = 1:Nnode
    indeces = (nn-1)*Nnode+1:nn*Nnode;
    e = [zeros((nn-1),1); 1; zeros(Nnode-nn,1)];
    Y = (e*e.')*Y_net;
    Yp(:,indeces) = .5.*(Y + Y');           % for injected active power
    Yq(:,indeces) = .5.*sqrt(-1).*(Y - Y'); % for injected reactive power
    M(:,indeces) = e*e.';                   % for voltage magnitude
end

Y = Y_net(2:end,2:end);
R = real(inv(Y));
X = imag(inv(Y));

Y_bar = Y_net(1,2:end);
Vnom = -inv(Y)*(Y_bar.').*MagPcc;

% for real part
g = real(R*diag(cos(angle(Vnom))./abs(Vnom)) - X*diag(sin(angle(Vnom))./Vnom));
h = real(X*diag(cos(angle(Vnom))./abs(Vnom)) + R*diag(sin(angle(Vnom))./Vnom));

for ss = 1:NumS
    for tt = tot_time:tot_time + horizon + 1
        PV_Errors_MPC_S(tt,:,ss) = zeros(7,1);
    end
end

A = cell(NumS,1);
for i = 1:NumS
    A{i} = zeros(tot_time+horizon,horizon,Nnode);
end

% prepare and check solar data forcast errors for distributionally robust
% optimization
for s = 1 : NumS
    for tt = 1:tot_time
        for th = 1 : horizon
            PP_pv_error = zeros(horizon,NumS);
            PP_pv_error = PV_Errors_MPC_S(tt,1:horizon+1,:); % within finite time horizon [2:horizon] have non-zero forecast error, 1 column is the real time, forecast error is zero.
            PP_pv_error = squeeze(PP_pv_error);
            PP_pv_error = PP_pv_error'; % P_pv (NumS,th)
            
            PV = PP_pv_error(s,th) + PV_Real24(tt + th);
            if PV < 0 % test if deterministic forecast + forecast errors are less than zero
                PV = 0;
            end
            for i = 1:Ninverters
                PPPV(i) = PV;
                if(PPPV(i) > maximum_S(i)./Sbase*1000)
                    PPPV(i) = maximum_S(i)./Sbase*1000;
                end
            end
            A{s}(tt,th,:) = zeros(Nnode,1);
            A{s}(tt,th,Node_PV) = PPPV;
        end
    end
end

%% parameters check before start
disp(['program and system setup done !']);
disp(['Please check the following settings before start optimization']);
disp(['Time Interval:   ', num2str(Timeinterval), ' mins'])
disp(['Finite time horizon: ', num2str(horizon)])
disp(['number of nodes:  ', num2str(Nnode)]);
disp(['number of energy storage devices:   ', num2str(Nbatt)]);
disp(['number of PV generators:  ', num2str(Ninverters)]);
disp(['number of scenarios for distributionally robust stochastic OPF:  ', num2str(NumS)]);
disp(['number of scenarios for Monte Carlo simulation:  ', num2str(NumMC)]);
disp(['range of weight factors:', 'start from ', num2str(min(rhoMatrix)),', up to ', num2str(max(rhoMatrix))]);
disp(['voltage lower bound: ', num2str(Vmin),', voltage upper bound: ', num2str(Vmax)]);


disp(['..................................................................................'])

disp(['press any keys to continue .......']);
pause


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  data-based distributionally robust stochastic OPF   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nrho = 1;
for i = 1:length(rhoMatrix)
    disp(['..................................................................................'])
    rho = rhoMatrix(i);
    curr_time = start_time; % current timestep
    beg_time = curr_time; % Store first time step
    
    %% Outputs %%%%%%%%%
    % Initialize
    V_f = zeros(Nnode, tot_time); % voltage solved by linear approximation
    V_a = zeros(Nnode, tot_time); % voltage solved by SDP relaxation
    B_f = zeros(Nbatt, tot_time); % batter state of charge
    Pb_f = zeros(Nbatt, tot_time); % battery charging power
    a_f = zeros(Ninverters, tot_time); % curtailment rate
    Q_f = zeros(Ninverters, tot_time); % reactive power
    Psub = zeros(1, tot_time);% power drawn from substations
    
    %% MPC
    while(curr_time <= tot_time) % Entire simulation time
        t1 = clock;
        cvx_solver MOSEK
        % cvx_solver SDPT3
        disp(['Current time: ' num2str(curr_time)]);
        disp(['Wasserestein distance is ', num2str(varepsilon)])
        disp(['Weighted factor is: ', num2str(rho)])
        
        cvx_begin quiet
        
        display(['solver for distritutionally robust stochastic OPF is  ', cvx_solver])
        
        variable a(Nnode-1, horizon) % Curtailed power
        variable PB(Nnode, horizon) % Dis/Charging from battery
        variable B(Nnode, horizon) % Energy level of battery
        variable Q(Nnode-1, horizon) % reactive power provided by inverter
        variable x(Ninverters, horizon) % Auxillary optimization variable
        
        if(curr_time > (60*6/Timeinterval) && curr_time < (60*18/Timeinterval))%
            variable sDROmax(Nnode-1, horizon, numS) % distributionally robust optimization auxillary varibles
            variable lambdamax(Nnode-1, horizon) % distributionally robust optimization auxillary varibles
            variable z(Nnode-1, horizon) % Auxillary optimization variable
            variable y(Nnode-1, horizon) % Auxillary optimization variable
        end
        
        if(curr_time >= (60*18/Timeinterval))
            variable y(Nnode-1, horizon) % Auxillary optimization variable
        end
        
        F = 0;
        for th = 1 : horizon
            d_l2 = dl{curr_time,th};
            
            F2 = 0;
            for s = 1 : NumS
                PP_pv = squeeze(A{s}(curr_time,th,:));
                F2 = F2 + f_i(curr_time+th-1)*sum((diag(a(:,th)))*(PP_pv(2:end))) + c_i(curr_time+th-1)*sum(max(0, P_l(2:end,curr_time+th-1) + PB(2:end,th) - (eye(Nnode-1)-diag(a(:,th)))*(PP_pv(2:end)))) + d_i(curr_time+th-1)*sum((max(0, -P_l(2:end,curr_time+th-1) - PB(2:end,th) + (eye(Nnode-1)-diag(a(:,th)))*(PP_pv(2:end)))));
            end
            
            % Reactive power cost
            F = F +  F2./NumS + e_i(curr_time+th-1).*sum(abs(Q(:,th)));
        end
        
        F3 = 0;
        if(curr_time > (60*6/Timeinterval) && curr_time < (60*18/Timeinterval))
            F3  = sum(sDROmax(:))/NumS + sum(lambdamax(:))*varepsilon;
        end
        
        F = F + F3;
        
        minimize  F
        subject to
        for t = curr_time : curr_time+horizon-1 % Optimization horizon
            % MPC constraints:
            % Get timestep within this horizon
            th = t-curr_time+1; % finite time horizon (MPC)
            
            disp(['HORIZON STEP: ' num2str(th)]);
            t3 = clock;
            if(th <= 6)
                QL(:,th) = Q_l(:,t);
                PL(:,th) = P_l(:,t);
            end
            
            %% Curtailed power limits
            a>=0;
            a<=1;
            
            % Set buses without PV to zero
            Q(Node_noPV-1,th) == 0;
            a(Node_noPV-1,th) == 0;
            PB(Node_nobatt,th) == 0;
            B(Node_nobatt,th) == 0;
            
            d_l2 = dl{curr_time,th};
            q_l = 0.*QL(2:end,th);
            p_l = 0.*PL(2:end,th);
            
            % Voltage constraints 6AM to 6PM
            p_ss = zeros(Nnode-1,1);
            
            if (curr_time >= 200)
                
                % loads
                q_l(:,th) = QL(2:end,th);
                p_l(:,th) = PL(2:end,th);
                
                summedG2 = zeros(length(Nnode-1), 1);
                summedG3 = zeros(Ninverters, 1);
                
                
                for s = 1:NumS
                    p_s = squeeze(A{s}(curr_time,th,2:end));
                    p_ss = p_ss + squeeze(A{s}(curr_time,th,2:end));
                    summedG3 = summedG3 + max(0, x(:,th) + ((1-a(Node_PV-1,th)).*p_s(Node_PV-1)).^2 + Q(Node_PV-1,th).^2 - p_s(Node_PV-1).^2);
                    gpn = abs(Vnom) + g*((eye(Nnode-1)-diag(a(:,th)))*p_s -p_l(:,th)-PB(2:end,th)) + h*(-q_l(:,th) + Q(:,th));
                end
                
                p_ss = p_ss./NumS;  % averaging solar forecast
                
                summedG2/NumS <= y(:,th)*epsilon; % voltage lower bound contraints (chance-constrained)
                summedG3/NumS <= x(:,th)*nu; % apparent power constraint of inverters (chance-constrained)
                
                % voltage upper bound constraints
                gpn = abs(Vnom) + g*((eye(Nnode-1)-diag(a(:,th)))*p_ss -p_l(:,th)-PB(2:end,th)) + h*(-q_l(:,th) + Q(:,th));
                gpn <= Vmax;
                
                % power factor constraint (sample average)
                abs(Q(Node_PV-1,th)) <= (eye(Ninverters)-diag(a(Node_PV-1,th)))*p_ss(Node_PV-1).*0.44; % Power Factor Constraints PF > 0.90
                
            end
            
            %% Upper voltage constraints
            % upper voltage constraints are only binding by distributionally
            % robust optimization during solar peak time
            if(curr_time > (60*6/Timeinterval) && curr_time < (60*18/Timeinterval))%
                
                summedG2 = zeros(length(Nnode-1), 1);
                summedG3 = zeros(Ninverters, 1);
                
                for s = 1 : NumS % Sum up each sample
                    q_l(:,th) = QL(2:end,th);
                    p_l(:,th) = PL(2:end,th);
                    
                    p_s = squeeze(A{s}(curr_time,th,2:end));
                    p_ss = p_ss + squeeze(A{s}(curr_time,th,2:end));
                    
                    gpn = abs(Vnom) + g*((eye(Nnode-1)-diag(a(:,th)))*p_s -p_l(:,th)-PB(2:end,th)) + h*(-q_l(:,th) + Q(:,th));
                    
                    summedG2 = summedG2 + max(0, y(:,th)  + Vmin - gpn); % voltage lower bound contraints
                    summedG3 = summedG3 + max(0, x(:,th) + ((1-a(Node_PV-1,th)).*p_s(Node_PV-1)).^2 + Q(Node_PV-1,th).^2 - p_s(Node_PV-1).^2); % apparent power constraint of inverters
                    
                    % distributionally robust optimization contraints
                    rho*(z(:,th) + gpn - z(:,th)*epsilon - Vmax) <= sDROmax(:,th,s);
                    rho*(-z(:,th)*epsilon) <=  sDROmax(:,th,s);
                end
                
                p_ss = p_ss./NumS;%(averaging solar forecast)
                
                summedG2/NumS <= y(:,th)*epsilon; % voltage lower bound contraints (chance-constrained)
                summedG3/NumS <= x(:,th)*nu; % apparent power constraint of inverters (chance-constrained)
                
                % power factor constraint (sample average)
                abs(Q(Node_PV-1,th)) <= (eye(Ninverters)-diag(a(Node_PV-1,th)))*p_ss(Node_PV-1).*0.44; % Power Factor Constraints PF > 0.95
                
                % distributionally robust optimization contraints
                for mm = 1:Nnode-1
                    norm([rho*g(mm,:).*(a(:,th)-1)', rho*g(mm,:).*PB(2:end,th)', -rho*h(mm,:).*Q(:,th)'], inf) <= lambdamax(mm,th);
                end
                norm(0, inf) <= lambdamax(:,th);
            end
            
            %% Battery dynamic and constraints
            % Link battery to previous actual timestep
            if(curr_time > beg_time && t == curr_time) % Beg of horizon
                B(Node_batt, 1) == Prev_B;
            elseif(curr_time == beg_time)
                B(Node_batt, 1) == B_0; % Initial SOC
            end
            if (th > 1)
                Bmin <= B(Node_batt,th) <= Bmax
            end
            
            if(th < horizon) % Link battery within horizon
                % Battery dynamics
                B(Node_batt,th+1) == B(Node_batt,th) + PB(Node_batt,th)*(Timeinterval/60);
                
                % Charging / Discharging Limits
                -ndis.*(B(Node_batt,th) - Bmin) <= PB(Node_batt,th)*(Timeinterval/60) <= nch.*(Bmax-B(Node_batt,th));
                PB(Node_batt,th) <= Bmax./10;
            end
            
            if(horizon == 1) % No MPC
                PB == 0;
                B == 0;
            end
        end
        
        t4 = clock;
        disp(['Problem formulation time :, ', num2str(etime(t4,t3))])
        cvx_end
        t5 = clock;
        disp(['Problem solving time : ', num2str(etime(t5,t1))])
        
        clear lhs
        disp(['CVX_STATUS: ', cvx_status]);
        if(isnan(cvx_optval) || strcmp(cvx_status, 'Infeasible'))
            disp('******** CVX FAILED TO CONVERGE *******')
        end
        Cost_Data(curr_time) = 0;
        CVaR_Data(curr_time) = 0;
        
        if(curr_time > (60*6/Timeinterval) && curr_time < (60*18/Timeinterval))%
            Cost_Data(curr_time) = cvx_optval - (sum(sDROmax(:))/NumS + sum(lambdamax(:))*varepsilon);
            CVaR_Data(curr_time) = (sum(sDROmax(:))/NumS + sum(lambdamax(:))*varepsilon)/rho;
            tot_cvx_cost(curr_time) = cvx_optval - (sum(sDROmax(:))/NumS + sum(lambdamax(:))*varepsilon);
        else
            Cost_Data(curr_time) = cvx_optval;
            tot_cvx_cost(curr_time) = cvx_optval;
        end
        tot_cvx_cputime(curr_time) = cvx_cputime;
        
        if(horizon > 1)
            Prev_B = B(Node_batt,2);
        else
            Prev_B = 0;
        end
        a2 = [zeros(1,size(a,2)); a]; % add one more column, which represent the a at the node 0, which does not have inverter around.
        
        % Store variables
        gpn = abs(Vnom) + g*((eye(Nnode-1)-diag(a(:,1)))*PV_Aggregate(2:end,curr_time) -PL(2:end,1)-PB(2:end,1)) + h*(-QL(2:end,1) + Q(:,1));
        
        V_f(:,curr_time) = [MagPcc gpn(:,1).'].';
        B_f(:,curr_time) = B(Node_batt,1);
        Pb_f(:, curr_time) = PB(Node_batt,1);
        a_f(:, curr_time) = a2(Node_PV,1);
        Q_f(:, curr_time) = Q(Node_PV-1,1);
        
        cvx_solver SDPT3 % Get actual voltages by SDP relaxation
        display(['The solver for actual voltages is  ', cvx_solver])
        [V_a(:,curr_time), Vsdp, success,cvx_status_VSDP] = getVSDP(Yp, Yq, P_l(:,curr_time),Q_l(:,curr_time),...
            PV_SDP(:,curr_time), MagPcc, ...
            Nnode, Node_PV, Node_noPV, Node_batt,...
            Ninverters, a_f(:,curr_time), Pb_f(:,curr_time),Q_f(:,curr_time));
        
        % Power drawn from substation
        Psub(curr_time) = real(trace(Yp(:,1:Nnode)*Vsdp));
        disp(['The CVX Status of get_VSDP:', cvx_status_VSDP])
        
        if(success == 0)
            disp('********* SDP did not satisfy rank 1 constraint **********');
        end
        disp('********* Start Monte Carlo Simultion, It might be a well ~~~~ **********');
        
        
        %% Monte Carlo Simluation
        for ss = 1:NumMC
            [V_a(:,curr_time), Vsdp, success,cvx_status_VSDP] = getVSDP(Yp, Yq, P_l(:,curr_time),Q_l(:,curr_time),...
                PV_MC_S(ss,:,curr_time), MagPcc, ...
                Nnode, Node_PV, Node_noPV, Node_batt,...
                Ninverters, a_f(:,curr_time), Pb_f(:,curr_time),Q_f(:,curr_time));
            
            % Power drawn from substation (Monte Carlo Simulation)
            Psub_MC(ss,curr_time) = real(trace(Yp(:,1:Nnode)*Vsdp));
            if(success == 0)
                disp('********* SDP did not satisfy rank 1 constraint **********');
            end
            
            % save the results from Monte Carlo simulation
            V_a_MC(ss,:,curr_time) = V_a(:,curr_time);
            success_MC(ss,curr_time) = success;
            cvx_status_VSDP_MC(ss,:,curr_time) = cvx_status_VSDP;
            
        end
        
        if(success_MC(:,curr_time) == 1 )
            disp('********* The Monte Carlo Simulation is Done at **********');
            display(['The weight factor is  ',  num2str(rho), ' , The Epsilon is  ', num2str(varepsilon),'  The current time is  ', num2str(curr_time)])
            disp('**********************************************************')
        end
        
        
        Data_Collection_CVaR(curr_time,Nrho) = CVaR_Data(curr_time);
        Data_Collection_Cost(curr_time,Nrho) = Cost_Data(curr_time);
        
        curr_time = curr_time+1; % Update the present time
        
        
    end
    MonteCarlo_Va(Nrho,:,:,:) = V_a_MC;
    MonteCarlo_success(Nrho,:,:) = success_MC;
    MonteCarlo_Psub(Nrho,:,:) = Psub_MC;
    
    Pb_tot=sum(Pb_f,1);
    
    a_tot = sum(a_f,1);
    totCurt = totPV(1:length(a_tot)).*a_tot/Ninverters;
    
    % save optimization results for selected weight factor \rho
    Actual_a_f(Nrho,:,:) = a_f;
    Actual_Pb_f(Nrho,:,:) = Pb_f;
    Actual_Voltage(Nrho,:,:) = V_a;
    Actual_Bf(Nrho,:,:) = B_f;
    Actual_Pb_tot(Nrho,:,:) = Pb_tot;
    Actual_Psub(Nrho,:,:) = Psub;
    Actual_totCurt(Nrho,:,:) = totCurt;
    Actual_Qf(Nrho,:,:) = Q_f';
    Actual_VoltageLinear(Nrho,:,:) = V_f; % Got the linearized voltage instead of SDP
    
    Nrho = Nrho + 1;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  results plot section   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1-st plot
figure;% CVaR value of voltage constraint violation under various weight factor
plot(Data_Collection_CVaR,'LineWidth',1.5)
legend('\rho = 1', '\rho = 5','\rho = 10', '\rho = 15','\rho = 20', '\rho = 30','\rho = 35', '\rho = 40','\rho = 100', '\rho = 200')
xlabel('Time Step')
ylabel('CVaR value of Voltage Constraints violations')
title(['CVaR value of Voltage Contraints ', char(949), ' = 0.000'])
% set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})


figure;% operational cost under various weight factor
plot(Data_Collection_Cost*Sbase/12/1000/100,'LineWidth',1.5)
legend('\rho = 1', '\rho = 5','\rho = 10', '\rho = 15','\rho = 20', '\rho = 30','\rho = 35', '\rho = 40','\rho = 100', '\rho = 200')
xlabel('Time Step')
ylabel('Operation Cost $/kWh')
title(['Operational Cost ' char(949) ' = 0.000'])
% set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})


%% 2-nd plot
% This plot visualize the voltage profiles solved by SDP relaxation, given
% all setpoints from distributionally robust stochastic OPF

figure;
for j = 1:length(rhoMatrix)
    subplot(floor(length(rhoMatrix)/2),2,j)
    plot(squeeze(Actual_Voltage(j,:,:)).') % \rho = 1
    title(['Actual Voltage \rho = ', num2str(rhoMatrix(j))]);
    xlabel('Time Step')
    ylabel('Voltage Magnitude in p.u.')
    % set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})
end

%% 3-rd plot
figure;
for j = 1:length(rhoMatrix)
    subplot(floor(length(rhoMatrix)/2),2,j)
    plot(squeeze(Actual_VoltageLinear(j,:,:))') % Actual voltages from SDP
    hold on
    plot(1:289,1.05.*ones(289),'k--','LineWidth',1,'color',[0.5 0.5 0.5])
    plot(1:289,0.95.*ones(289),'k--','LineWidth',1,'color',[0.5 0.5 0.5])
    title(['Linearized Voltage, \rho = ', num2str(rhoMatrix(j)),' , ', char(949),  '  = ', num2str(varepsilon)]);
    xlabel('Time Step')
    ylabel('p.u.')
    % set(gca,'xtick',[1 36 72 108 144 180 204 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','5 PM','6 PM','9 PM','12 AM'})
    set(gca,'FontSize',11)
    axis([72 204 0.93 1.07])
    grid
end


%% 4-th plot
% This plot visualizes the power drawn from substation under various weight
% factor \rho
figure;
for j = 1:length(rhoMatrix)
    plot(squeeze(Actual_Psub(j,:,:)).*Sbase/1000, 'lineWidth', 2) % rho = 1
    hold on
end
xlabel('Time Step')
ylabel('kW')
title(['Power Drawn from Substation  ', char(949), '  = 0.0000'])
legend('\rho = 1', '\rho = 15','\rho = 20', '\rho = 30','\rho = 35','\rho = 100','\rho = 200')
% set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})
set(gca,'FontSize',13)


%% 5-th plot
% This plot visualize the total solar power curtailment under various
% weight factor \rho
figure;
for j = 1:length(rhoMatrix)
    plot(squeeze(Actual_totCurt(j,:,:).*Sbase/1000), 'lineWidth', 2)
    hold on
end
xlabel('Time Step')
ylabel('kW')
title(['Total Curtailed Power  ', char(949), '  = ', num2str(varepsilon)])
legend('\rho = 1', '\rho = 15','\rho = 20', '\rho = 30','\rho = 35','\rho = 100','\rho = 200')
% set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})
set(gca,'FontSize',13)

%% 6-th plot
% This plot visualize the total solar power injection under various weight
% factor \rho
figure;
for j = 1:length(rhoMatrix)
    plot(((totPV(1:length(a_tot))' - squeeze(Actual_totCurt(j,:,:))).*Sbase/1000), 'lineWidth', 2)
    hold on
end
xlabel('Time Step')
ylabel('kW')
title(['Total PV Injection  ', char(949), '  = ', num2str(varepsilon)])
legend('\rho = 1', '\rho = 15','\rho = 20', '\rho = 30','\rho = 35','\rho = 100','\rho = 200')
% set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})
set(gca,'FontSize',13)

%% 7-th plot
% This plot visualize the battery charging power under various weight
% factor \rho
figure;
for j = 1:length(rhoMatrix)
    plot(squeeze(Actual_Pb_tot(1,:,:).*Sbase/1000), 'lineWidth', 2)
    hold on
end
title(['Power Dis/charged from Storage ', char(949), '  = ', num2str(varepsilon)]);
xlabel('Time Step')
ylabel('kW.')
legend('\rho = 1', '\rho = 5','\rho = 10', '\rho = 15','\rho = 20', '\rho = 30','\rho = 35', '\rho = 40','\rho = 100', '\rho = 200')
%     set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})

%% 8-th plot
% This plot visualize the reactive power compensation under various weight
% factor \rho
figure;
for j = 1:length(rhoMatrix)
    subplot(floor(length(rhoMatrix)/2),2,j)
    plot(squeeze(Actual_Qf(1,:,:).*Sbase/1000), 'lineWidth', 2)
    title(['Reactive Power Compensation in kVar \rho = ', num2str(rhoMatrix(j))]);
    % xlabel('Time Step')
    set(gca,'xtick',[1 36 72 108 144 180 216 252 288],'xticklabel',{'12 AM','3 AM','6 AM','9 AM','12 PM','3 PM','6 PM','9 PM','12 AM'})
    ylabel('kvar')
end

disp('Finished!')

toc
