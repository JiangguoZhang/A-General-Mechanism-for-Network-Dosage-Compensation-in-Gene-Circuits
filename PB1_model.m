function gal = PB1_model(v,z, functype)
    % Parameters
    theta_2 = 1500; %1500;
    theta_3 = 1500; %1500;
    theta_4 = 100; %100;
    theta_80 = 1500; %1500;
    lambda = 0.2;
    gamma = 0.46;

    % Calculate functional
    if functype == 1
        F = functional_form_1(z(1), z(2), z(3), z(4)); %normal
    end
    if functype == 2
        F = functional_form_2(z(1), z(2), z(3), z(4)); %engineered
    end
        
    %ODES FOR gal2, gal3, gal4, gal80
    gal = zeros(4,1);
    gal(1) = theta_2 * (lambda + (1 - lambda) *  F) - gamma * z(1); %gal2
    gal(2) = theta_3 * (lambda + (1 - lambda) *  F) - gamma * z(2); %gal3
    gal(3) = theta_4  - gamma * z(3); %gal4
    gal(4) = theta_80 * (lambda + (1 - lambda) *  F) - gamma * z(4); %gal80
end