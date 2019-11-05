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
beta = 5;
s_a = 100;
s_i = 10;
theta_a = 1500 ; %constant
theta_i = 1500 ; %constant

%% Find the changing values of alpha for the 2 gene system
concentrations=[];
for conc=logspace(-2,0,7)
    g = conc;
    theta_a = 1500;
    theta_i = 1500;
    [v1 z1] = ode45(@myode,[0 12],[1 1]);
    theta_a = 1500;
    theta_i = 750;
    [v2 z2] = ode45(@myode,[12 12.5],[z1(end,1) z1(end,2)]);
    theta_a = 1500;
    theta_i = 1500;
    [v3 z3] = ode45(@myode,[12.5 24],[z2(end,1) z2(end,2)]);
    v=[v1;v2;v3];
    z=[z1;z2;z3];
    fraction=( 1 ./ ( 1 + ( s_i .* (z(:,2)) ./ (1 + (s_a .* g .* z(:,1)).^alpha)).^beta ) );
    %subplot(1,2,1);
    plot(v,fraction);
    hold on;
    %subplot(1,2,2);
    %plot(v,z(:,2));
    hold on;
    concentrations = [concentrations; join(["g=",num2str(conc)],"")];
end

%semilogx(concentrations,cell_frac);
%hold on;

% subplot(1,2,1);
% title('activator');
xlabel('time/h');
ylabel("ON state fraction");
% ylabel('Concentration');
% legend(concentrations,'Location','best');
% subplot(1,2,2);
% title('inhibitor');
% xlabel('time/h');
% ylabel('Concentration');
legend(concentrations,'Location','best');



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