
function functional = functional_form_1(x_2t, x_3t, x_80t, x_4t)
    global gg;
    global alpha;
    global beta;
    global u; 
    
    %% Galactose count 
    gt = gg;
 
    %% Parameters (Chosen from Supp. info; page 11)
    n = 1.6;
    k_4 = 8; %2;
    k_80 = 8; 
    v = 1.3;
    u_m = u; %0.9;
    k_g = 0.052;
    k_2 = 1700; %600;
    k_3 = 4.1;
    alpha_m = alpha; % 0.85;
    beta_m = beta; %6.1;
    h =2.5;
  
    %% Equations for the functional form (as in the supplemental info)
    g = ( gt / (x_2t / k_2)^(-u_m) );
    x_80 = (x_80t /  (1 + (g / k_g)^(-v)) );
    x_3 = (x_3t / ( 1 + (x_80 / k_80)^beta_m ) );
    x_4 = x_4t / (1 + (x_3 / k_3)^(-alpha_m));
    rho = (x_4 / k_4);
    functional = (1 / (1 + (1/rho)) );
end
