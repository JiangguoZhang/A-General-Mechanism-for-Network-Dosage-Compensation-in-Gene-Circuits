clc; close all; clear;
% PBL model for "A General Mechanism for Network-Dosage Compensationin Gene
% Circuits"

% TEAM: Dreycey Albin, Jiangguo Zhang, Cole Grandel, Kiara Reyes Gamas

% Assign the global variables for the system
global g;
global alpha;
global beta;
global s_a;
global s_i;
global theta_a;
global theta_i;

alpha = 1;
beta = 1;
s_a = 100;
s_i = 10;
theta_a = 750 ; %constant
theta_i = 750 ; %constant

%% Find the changing values of alpha for the 2 gene system

alpha_values=[];
for alph=0.5:0.5:5
    alpha=alph;
    alpha_values = [alpha_values; join(["\alpha=",num2str(alph)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_i * (z(end,2)) / (1 + (s_a * g * z(end,1))^alpha))^beta ) );
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,1)
    semilogx(concentrations,cell_frac, 'color', [1-alph/5,0,alph/5]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing \alpha)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(alpha_values), 'location','eastoutside');

alpha=1;
beta_values = [];
for bet=0.5:0.5:5
    beta = bet;
    beta_values = [beta_values; join(["\beta=",num2str(beta)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_i * (z(end,2)) / (1 + (s_a * g * z(end,1))^alpha))^beta ) );
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,2)
    semilogx(concentrations,cell_frac, 'color', [1-bet/5,0,bet/5]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing \beta)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(beta_values), 'location','eastoutside');

beta = 1;
si_values = [];
for si=logspace(-3,1,10)
    s_i = si;
    si_values = [si_values; join(["s_i=",num2str(si)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_i * (z(end,2)) / (1 + (s_a * g * z(end,1))^alpha))^beta ) );
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,3)
    semilogx(concentrations,cell_frac, 'color', [1-(log10(si)+3)/4,0,(log10(si)+3)/4]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing s_i)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(si_values), 'location','eastoutside');

s_i = 10;
sa_values = [];
for sa=logspace(-2,2,10)
    s_a = sa;
    sa_values = [sa_values; join(["s_a=",num2str(sa)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_i * (z(end,2)) / (1 + (s_a * g * z(end,1))^alpha))^beta ) );
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,4)
    semilogx(concentrations,cell_frac, 'color', [1-(log10(sa)+2)/4,0,(log10(sa)+2)/4]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing s_a)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(sa_values), 'location','eastoutside');

%% Find the changing values of alpha for the 2 gene system
beta = 1;
s_a = 100;
s_i = 10;
alpha_values=[];
for alph=0.5:0.5:5
    alpha=alph;
    alpha_values = [alpha_values; join(["\alpha=",num2str(alph)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode_2,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_a * g * (z(end,2)) / (s_i * z(end,1))^beta)^(-alpha)));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,1)
    semilogx(concentrations,cell_frac, 'color', [1-alph/5,0,alph/5]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing \alpha)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(alpha_values), 'location','eastoutside');

alpha=1;
beta_values = [];
for bet=0.5:0.5:5
    beta = bet;
    beta_values = [beta_values; join(["\beta=",num2str(beta)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode_2,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_a * g * (z(end,2)) / (s_i * z(end,1))^beta)^(-alpha)));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,2)
    semilogx(concentrations,cell_frac, 'color', [1-bet/5,0,bet/5]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing \beta)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(beta_values), 'location','eastoutside');

beta = 1;
si_values = [];
for si=logspace(-3,1,10)
    s_i = si;
    si_values = [si_values; join(["s_i=",num2str(si)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode_2,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_a * g * (z(end,2)) / (s_i * z(end,1))^beta)^(-alpha)));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,3)
    semilogx(concentrations,cell_frac, 'color', [1-(log10(si)+3)/4,0,(log10(si)+3)/4]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing s_i)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(si_values), 'location','eastoutside');

s_i = 10;
sa_values = [];
for sa=logspace(-2,2,10)
    s_a = sa;
    sa_values = [sa_values; join(["s_a=",num2str(sa)],"")];
    counter = 1;
    steadyStates=[];
    concentrations=[];
    cell_frac=[];
    for conc=logspace(-2,0,100)
        g = conc;
        [v z] = ode45(@myode_2,[0 29],[1 1]);
        steadyStates = [steadyStates; z(end,1)];
        F = ( 1 / ( 1 + ( s_a * g * (z(end,2)) / (s_i * z(end,1))^beta)^(-alpha)));
        cell_frac = [cell_frac; F*100];
        concentrations = [concentrations; g];
        counter = counter +1;
    end
    subplot(2,2,4)
    semilogx(concentrations,cell_frac, 'color', [1-(log10(sa)+2)/4,0,(log10(sa)+2)/4]);
    hold on;
end
ylim([0, 100])
title('Response to galactose concentration (changing s_a)');
xlabel('concentration of galactose');
ylabel('Fraction of ON cells [%]');
legend(cellstr(sa_values), 'location','eastoutside');



function dz = myode(v,z)
global g;
global alpha;
global beta;
global s_a;
global s_i;
global theta_a;
global theta_i;
alpha_m = alpha;
beta_m = beta;
lambda = 0.2 ; %constant
gamma = 0.46 ; %constant
g_model = g;

dz = zeros(2,1); 

dz(1) = theta_a * ( lambda + (1 - lambda) * ( 1 / ( 1 + ( s_i * (z(2)) / (s_a * (z(1)) * g_model)^alpha_m)^beta_m ) )) - gamma * (z(1));
dz(2) = theta_i * ( lambda + (1 - lambda) * ( 1 / ( 1 + ( s_i * (z(2)) / (s_a * (z(1)) * g_model)^alpha_m)^beta_m ) )) - gamma * (z(2));

end

function dz = myode_2(v,z)
global g;
global alpha;
global beta;
global s_a;
global s_i;
global theta_a;
global theta_i;
alpha_m = alpha;
beta_m = beta;
lambda = 0.2 ; %constant
gamma = 0.46 ; %constant
g_model = g;

dz = zeros(2,1); 

dz(1) = theta_a * ( lambda + (1 - lambda) * ( 1 / ( 1 + ( s_i * (z(2) *g_model) / (s_a * (z(1)) )^beta_m)^(-alpha_m) ) )) - gamma * (z(1));
dz(2) = theta_i * ( lambda + (1 - lambda) * ( 1 / ( 1 + ( s_i * (z(2)* g_model) / (s_a * (z(1)) )^beta_m)^(-alpha_m) ) )) - gamma * (z(2));

end

%{
NOTES:

EQUATIONS FOR 2 GENE MODEL:
(Partition Function) Z = (1 + (((z(1))^2) / kx) +  (((z(2))^2) / ky) + (((z(2))^2) * (z(2))^2)/ (ky * kx))
F_1 = ( 1 / ( 1 + ( s_i * (z(2)) / (s_a * (z(1)) * g_model)^alpha_m)^beta_m ) )
F_2 = ( 1 / ( 1 + ( s_i * (z(2)) / (s_a * (z(1)) * g_model)^beta_m)^(-alpha_m) ) )
dz(1) = theta_a * ( lambda + (1 - lambda) * F) - gamma * (z(1));
dz(2) = theta_a * ( lambda + (1 - lambda) * F) - gamma * (z(2));
%}

