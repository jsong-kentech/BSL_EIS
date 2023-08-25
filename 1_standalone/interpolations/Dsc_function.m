function Dsc_est = Dsc_function( theta0c,T)
% This function calculates the solid-state diffusivity of anode active
% material. 
% Measured data is from 'p_am1_Ds_HC4088HC4101HC4109.txt'

%% Measured data
    data = [0.3	45287.92969	1.64077E-14
            0.31	44653.02566	1.66162E-14
            0.32	43819.77494	1.69333E-14
            0.33	42818.10435	1.7349E-14
            0.34	41677.94197	1.78536E-14
            0.35	40429.21511	1.84372E-14
            0.36	39101.8531	1.90887E-14
            0.37	37725.78091	1.97957E-14
            0.38	36330.92884	2.05433E-14
            0.39	34947.22293	2.13146E-14
            0.4	33604.59106	2.20894E-14
            0.41	32332.96251	2.28451E-14
            0.42	31162.26261	2.35557E-14
            0.43	30122.42095	2.41933E-14
            0.44	29242.19049	2.47283E-14
            0.45	28541.67412	2.51347E-14
            0.46	28037.10698	2.5389E-14
            0.47	27744.70461	2.54696E-14
            0.48	27680.68206	2.5358E-14
            0.49	27861.25611	2.50403E-14
            0.5	28302.64293	2.45077E-14
            0.51	29021.05839	2.37586E-14
            0.52	30002.19216	2.28119E-14
            0.53	31112.73282	2.17406E-14
            0.54	32172.94382	2.06196E-14
            0.55	32992.21614	1.95111E-14
            0.56	33379.79473	1.8466E-14
            0.57	33144.92854	1.75247E-14
            0.58	32096.86252	1.6719E-14
            0.59	30132.45703	1.60477E-14
            0.6	27806.53754	1.53142E-14
            0.61	25959.24421	1.42792E-14
            0.62	25072.13596	1.28733E-14
            0.63	25216.4873	1.11931E-14
            0.64	26230.60827	9.41022E-15
            0.65	27955.12679	7.70737E-15
            0.66	30077.24602	6.23031E-15
            0.67	32399.66187	5.02447E-15
            0.68	35210.45517	4.04813E-15
            0.69	38390.36841	3.28871E-15
            0.7	40795.11995	2.77351E-15
            0.71	41817.65944	2.47393E-15
            0.72	42339.52197	2.31058E-15
            0.73	43206.26765	2.22632E-15
            0.74	44530.46354	2.17654E-15
            0.75	46260.7889	2.12395E-15
            0.76	47992.10898	2.05434E-15
            0.77	49417.87855	1.97987E-15
            0.78	50446.53612	1.9128E-15
            0.79	51303.59877	1.85482E-15
            0.8	52421.0795	1.80527E-15
            0.81	54323.61663	1.7602E-15
            0.82	57412.55052	1.70235E-15
            0.83	61769.082	1.61148E-15
            0.84	67363.76582	1.47773E-15
            0.85	74359.51726	1.30885E-15
            0.86	83163.73189	1.11815E-15
            0.87	94133.19937	9.12E-16
            0.88	106655.8489	6.94574E-16
            0.89	119420.9053	4.86979E-16
            0.9	130091.6361	3.1927E-16
            0.91	136061.1173	1.73583E-16
            0.92	136137.6562	1.01391E-16
            0.93	129545.1365	6.00497E-17
            0.94	116031.8862	3.74034E-17
            0.95	95782.41185	2.54409E-17];
%% calculates
    D0 = interp1(data(:,1),data(:,3),theta0c,'linear','extrap'); % [m2/s] prefactor interpolation (D at T=298.15K)
    Ea = interp1(data(:,1),data(:,2),theta0c,'linear','extrap'); % [J/mol] Ea interpolation.
    R = 8.314; % [J.mol/K]
    Dsc_est = D0.*exp(-Ea/R*(1/T - 1/298.15));
    
end
