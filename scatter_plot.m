clc; close all; clear;
%%
syms lg;
gtyp = 1/10;%si/sa*ti/ta;
k = 5.5;
f=@(lg) 1./(1+(gtyp./exp(lg)).^k);

% x = logspace(-2,0,100);
% y = f(x);
% semilogx(x,y);
areaw = int(f,lg,-2,0);
areaw = eval(areaw);

%%
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
theta_a = 1500 ; %constant
theta_i = 1500 ; %constant

alpha=1;
beta=5;
s_i=10;
s_a=100;
theta_a = 1500;
theta_i = 1500;
areaw = 0;
steps=200;
xrange = 2;
for conc=logspace(-2,0,steps)
    g = conc;
    [v, z] = ode45(@myode,[0 24],[1 1]); %diploid
    fraction = ( 1 / ( 1 + ( s_i * (z(end,2)) / (1 + (s_a * g * z(end,1))^alpha))^beta ) );
    areaw = areaw + fraction / steps * xrange;
end

%% Find the changing values of alpha for the 2 gene system
steps=200;
xrange = 2;
plotdata = [];
for alph=0.5:0.5:5
    alpha=alph
    for bet=0.5:0.5:5
        beta=bet
        for si=logspace(-3,1,5)
            s_i=si;
            for sa=logspace(-2,2,5)
                s_a=sa;
                theta_a = 1500;
                theta_i = 1500;
                area2 = 0;
                for conc=logspace(-2,0,steps)
                    g = conc;
                    [v, z] = ode45(@myode,[0 24],[1 1]); %diploid
                    fraction = ( 1 / ( 1 + ( s_i * (z(end,2)) / (1 + (s_a * g * z(end,1))^alpha))^beta ) );
                    area2 = area2 + fraction / steps * xrange;
                end
                theta_a = 750;
                theta_i = 750;
                area1 = 0;
                for conc=logspace(-2,0,steps)
                    g = conc;
                    [v, z] = ode45(@myode,[0 24],[1 1]); %haploid
                    fraction = ( 1 / ( 1 + ( s_i * (z(end,2)) / (1 + (s_a * g * z(end,1))^alpha))^beta ) );
                    area1 = area1 + fraction / steps * xrange;
                end
                plotdata = [plotdata; [abs(area2-area1), abs(area2-areaw),alph,bet,si,sa]];
            end
        end
    end
end
%%
scatter(plotdata(:,2)*5, plotdata(:,1)*5,'*');
set(gca,'xscale','log');
set(gca,'yscale','log');
%xlim([0.01,10]);
%ylim([0.01,100]);
xlabel("penalty to wildtype inducibility [a.u.]");
ylabel("penalty to network dosage invariance [a.u.]");

%%
aha = [plotdata(:,1).*plotdata(:,2),plotdata];
newmat = aha(:,1)<0.01;
r1 = aha(:,2)<0.5;
r2 = aha(:,3)<0.5;
newmat = newmat .* r1 .* r2 == 1;
fin = aha(newmat,:);

%%
subplot(2,2,1);
histogram(fin(:,4));
title('\alpha');
subplot(2,2,2);
histogram(fin(:,5));
title('\beta');
subplot(2,2,3);
histogram(log10(fin(:,6)));
title('s_i');
subplot(2,2,4);
histogram(log10(fin(:,7)));
title('s_a');


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

