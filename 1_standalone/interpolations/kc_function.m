function [ kc ] = kc_function(theta0c,T,ctc)
% This function calculates the exchange current density of anode material
% data is from "p_am1_i0_preexp_Ea_22FKV08_HC2418-19.txt"
%% data [ x    i_0_preexp    Ea]
data = [0.3228	5.2931	43109.7569
        0.3290	5.9810	44878.9315
        0.3500	5.9030	44019.4454
        0.4000	5.5784	45105.8757
        0.4500	5.1843	45270.3738
        0.5000	4.8107	46118.0210
        0.5500	4.7017	46095.5915
        0.6000	4.8911	44541.9238
        0.6500	4.8922	44354.0326
        0.7000	4.6196	45358.9255
        0.7500	3.9369	48624.7728
        0.8000	2.6872	54390.1304
        0.8500	1.6266	63227.0378
        0.8809	0.8933	59213.0862
        0.8907	0.8126	55822.1118
        0.9220	0.7414	39083.2541];

%% Interpolation
    i0c_preexp = interp1(data(:,1),data(:,2),theta0c,'linear','extrap');
    Ea = interp1(data(:,1),data(:,3),theta0c,'linear','extrap');
    R = 8.314;
    F = 96487;
    kc = i0c_preexp*exp(-Ea/R*(1/T - 1/298.15))/F/ctc;
    
end

