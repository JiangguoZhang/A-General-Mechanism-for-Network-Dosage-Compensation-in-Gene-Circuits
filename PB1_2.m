clc; close all; clear all;
%PB1 model extension 
% TEAM: Dreycey Albin, Jiangguo Zhang, Cole Grandel, Kiara Reyes Gamas

%{
The following is a series of scripts for generating the different figures
for PB1. These scripts use the files and functions 'functional_form_1.m',
'functional_form_2.m', and 'PB1_model.m'. 
%}

% Assign the global variables for the system
global gg;
global alpha;
global beta;
global u;

%%%
% Looking at how GAL_i concentrations change per model
%%%
%plotting parameters
yhigh = 3500; %height of y axis
NF =50; % Noise Factor
% Keeping alpha and beta as is in supp. mat.
alpha= 0.85; %given in supp.
beta= 6.1; %given in supp.
u = 0.9;
gg=0.01;
% plotting the concentrations of the GAL proteins for the regular network
functype = 1;
subplot(1,2,1);
[v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
plot(v,z(:,1)+NF*randn(length(z)),'-o')
ylim([0, yhigh]);
hold on
plot(v,z(:,2)+NF*randn(length(z)),'-o')
ylim([0, yhigh]);
hold on
plot(v,z(:,3)+NF*randn(length(z)),'-o', 'Color', 'black')
ylim([0, yhigh]);
hold on
plot(v,z(:,4)+NF*randn(length(z)),'-o')
ylim([0, yhigh]);
title('Time solution NEG feedback');
xlabel('Time (min)');
ylabel('[GAL_i]');
%legend('GAL2','GAL3','GAL4','GAL80')
z(end,1)
z(end,2)
z(end,3)
z(end,4)

% plotting the concentrations of the GAL proteins for the modified network
functype = 2; 
subplot(1,2,2);
[v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
plot(v,z(:,1)+NF*randn(length(z)),'-o')
ylim([0, yhigh]);
hold on
plot(v,z(:,2)+NF*randn(length(z)),'-o')
ylim([0, yhigh]);
hold on
plot(v,z(:,3)+NF*randn(length(z)),'-o', 'Color', 'black')
ylim([0, yhigh]);
hold on
plot(v,z(:,4)+NF*randn(length(z)),'-o')
ylim([0, yhigh]);
title('Time solution POS feedback');
xlabel('Time (min)');
ylabel('[GAL_i]');
legend('GAL2','GAL3','GAL4','GAL80')
figure();



%%%
% Measuring effects of parameters
%%%
% Graph showing the decrease in concentration wrt galactose function #2
% Graph showing the decrease in concentration wrt galactose function #2
for bet= 0:0.1:10 %0.8:0.02:1.2
    alpha=bet;
    beta = 1;
    u = 0.9;
    cell_frac=[];
    steadyStates=[];
    concentrations=[];
    counter=0;
    for conc=0:0.01:1
        gg = conc;
        functype=1;
        [v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
        steadyStates = [steadyStates; z(end,3)];
        F = functional_form_1(z(end,1), z(end,2), z(end,3), z(end,4));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; gg];
        counter = counter + 1;
    end
    subplot(2,3,1);
    semilogx(concentrations,cell_frac,'-o'); %, 'Color', 'black');
    hold on
end
title('Response to galactose concentration (NEG FEEDBACK TOPOLOGY) varying alpha');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');

% Graph showing the decrease in concentration wrt galactose function #2
for bet= 0:0.1:10 %0.8:0.02:1.2
    beta=bet;
    alpha = 0.85;
    u = 0.9;
    cell_frac=[];
    steadyStates=[];
    concentrations=[];
    counter=0;
    for conc=0:0.01:1
        gg = conc;
        functype=1;
        [v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
        steadyStates = [steadyStates; z(end,3)];
        F = functional_form_1(z(end,1), z(end,2), z(end,3), z(end,4));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; gg];
        counter = counter + 1;
    end
    subplot(2,3,2);
    semilogx(concentrations,cell_frac,'-o'); %, 'Color', 'black');
    hold on
end
title('Response to galactose concentration (NEG FEEDBACK TOPOLOGY) varying beta');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');

% Graph showing the decrease in concentration wrt galactose function #2
for bet= 0:0.1:10 %0.8:0.02:1.2
    beta=6;
    alpha = 0.85;
    u = bet;
    cell_frac=[];
    steadyStates=[];
    concentrations=[];
    counter=0;
    for conc=0:0.01:1
        gg = conc;
        conc
        functype=1;
        [v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
        steadyStates = [steadyStates; z(end,3)];
        F = functional_form_1(z(end,1), z(end,2), z(end,3), z(end,4));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; gg];
        counter = counter + 1;
    end
    subplot(2,3,3);
    semilogx(concentrations,cell_frac,'-o'); %, 'Color', 'black');
    hold on
end
title('Response to galactose concentration (NEG FEEDBACK TOPOLOGY) varying u');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');

% Graph showing the decrease in concentration wrt galactose function #2
for bet= 0:0.1:10 %0.8:0.02:1.2
    alpha=bet;
    beta = 6;
    u = 0.9;
    cell_frac=[];
    steadyStates=[];
    concentrations=[];
    counter=0;
    for conc=0:0.01:1
        gg = conc;
        conc
        functype=2;
        [v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
        steadyStates = [steadyStates; z(end,3)];
        F = functional_form_2(z(end,1), z(end,2), z(end,3), z(end,4));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; gg];
        counter = counter + 1;
    end
    subplot(2,3,4);
    semilogx(concentrations,cell_frac,'-o'); %, 'Color', 'black');
    hold on
end
title('Response to galactose concentration (POS FEEDBACK TOPOLOGY) varying alpha');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');

% Graph showing the decrease in concentration wrt galactose function #2
for bet= 0:0.1:10 %0.8:0.02:1.2
    beta=bet;
    alpha = 0.85;
    u = 0.9;
    cell_frac=[];
    steadyStates=[];
    concentrations=[];
    counter=0;
    for conc=0:0.01:1
        gg = conc;
        conc
        functype=2;
        [v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
        steadyStates = [steadyStates; z(end,3)];
        F = functional_form_2(z(end,1), z(end,2), z(end,3), z(end,4));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; gg];
        counter = counter + 1;
    end
    subplot(2,3,5);
    semilogx(concentrations,cell_frac,'-o'); %, 'Color', 'black');
    hold on
end
title('Response to galactose concentration (POS FEEDBACK TOPOLOGY) varying beta');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');

% Graph showing the decrease in concentration wrt galactose function #2
for bet= 0:0.1:10 %0.8:0.02:1.2
    beta=6;
    alpha = 0.85;
    u = bet;
    cell_frac=[];
    steadyStates=[];
    concentrations=[];
    counter=0;
    for conc=0:0.01:1
        gg = conc;
        conc
        functype=2;
        [v z] = ode45(@(v,z) PB1_model(v,z,functype),[0 24],[1 1 1 1]);
        steadyStates = [steadyStates; z(end,3)];
        F = functional_form_2(z(end,1), z(end,2), z(end,3), z(end,4));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; gg];
        counter = counter + 1;
    end
    subplot(2,3,6);
    semilogx(concentrations,cell_frac,'-o'); %, 'Color', 'black');
    hold on
end
title('Response to galactose concentration (POS FEEDBACK TOPOLOGY) varying u');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');



%%%%
% 
%%%%