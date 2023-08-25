function [De_est] = De_function( c_e, T )
% Interpolates the diffusivity of electrolytes from a measured data
% Measured data is from 'De_comsol.txt'

%% data
c_e_sample = [0.0250000000000000,0.0500000000000000,0.100000000000000,0.150000000000000,0.200000000000000,0.250000000000000,0.300000000000000,0.350000000000000,0.400000000000000,0.450000000000000,0.500000000000000,0.550000000000000,0.600000000000000,0.650000000000000,0.700000000000000,0.750000000000000,0.800000000000000,0.850000000000000,0.900000000000000,0.950000000000000,1,1.05000000000000,1.10000000000000,1.15000000000000,1.20000000000000,1.25000000000000,1.30000000000000,1.35000000000000,1.40000000000000,1.45000000000000,1.50000000000000,1.55000000000000,1.60000000000000,1.65000000000000,1.70000000000000,1.75000000000000,2,2.50000000000000,2.75000000000000,3];
T_sample = [243  253	263	273	283	293	298	303	313	323];
De_sample = [5.1800e-13	2.0700e-11	9.4900e-11	2.1800e-10	3.6700e-10	5.2600e-10	6.0600e-10	6.8400e-10	8.3500e-10	9.7800e-10
4.7200e-13	1.9900e-11	9.2500e-11	2.1300e-10	3.6100e-10	5.1800e-10	5.9600e-10	6.7300e-10	8.2300e-10	9.6400e-10
3.9100e-13	1.8400e-11	8.7700e-11	2.0500e-10	3.4800e-10	5.0100e-10	5.7700e-10	6.5300e-10	7.9900e-10	9.3600e-10
3.2100e-13	1.7000e-11	8.3200e-11	1.9600e-10	3.3600e-10	4.8400e-10	5.5900e-10	6.3300e-10	7.7500e-10	9.1000e-10
2.6300e-13	1.5600e-11	7.8900e-11	1.8800e-10	3.2400e-10	4.6900e-10	5.4200e-10	6.1300e-10	7.5300e-10	8.8400e-10
2.1300e-13	1.4400e-11	7.4800e-11	1.8000e-10	3.1200e-10	4.5300e-10	5.2400e-10	5.9500e-10	7.3000e-10	8.5800e-10
1.7200e-13	1.3200e-11	7.0800e-11	1.7300e-10	3.0100e-10	4.3900e-10	5.0800e-10	5.7600e-10	7.0900e-10	8.3400e-10
1.3700e-13	1.2100e-11	6.7000e-11	1.6600e-10	2.9000e-10	4.2400e-10	4.9200e-10	5.5900e-10	6.8800e-10	8.1000e-10
1.0900e-13	1.1100e-11	6.3400e-11	1.5900e-10	2.8000e-10	4.1000e-10	4.7600e-10	5.4100e-10	6.6800e-10	7.8700e-10
8.5700e-14	1.0100e-11	6.0000e-11	1.5200e-10	2.6900e-10	3.9700e-10	4.6100e-10	5.2500e-10	6.4800e-10	7.6500e-10
6.6800e-14	9.2400e-12	5.6700e-11	1.4600e-10	2.6000e-10	3.8400e-10	4.4600e-10	5.0900e-10	6.2900e-10	7.4300e-10
5.1500e-14	8.4300e-12	5.3600e-11	1.4000e-10	2.5000e-10	3.7100e-10	4.3200e-10	4.9300e-10	6.1000e-10	7.2100e-10
3.9400e-14	7.6700e-12	5.0600e-11	1.3400e-10	2.4100e-10	3.5900e-10	4.1800e-10	4.7800e-10	5.9200e-10	7.0100e-10
2.9700e-14	6.9700e-12	4.7800e-11	1.2800e-10	2.3200e-10	3.4700e-10	4.0500e-10	4.6300e-10	5.7500e-10	6.8100e-10
2.2200e-14	6.3200e-12	4.5100e-11	1.2200e-10	2.2400e-10	3.3500e-10	3.9200e-10	4.4800e-10	5.5800e-10	6.6100e-10
1.6300e-14	5.7300e-12	4.2500e-11	1.1700e-10	2.1600e-10	3.2400e-10	3.8000e-10	4.3400e-10	5.4100e-10	6.4200e-10
1.1900e-14	5.1800e-12	4.0100e-11	1.1200e-10	2.0800e-10	3.1400e-10	3.6700e-10	4.2100e-10	5.2500e-10	6.2400e-10
8.4800e-15	4.6700e-12	3.7800e-11	1.0700e-10	2.0000e-10	3.0300e-10	3.5600e-10	4.0800e-10	5.1000e-10	6.0600e-10
5.9700e-15	4.2100e-12	3.5500e-11	1.0200e-10	1.9200e-10	2.9300e-10	3.4400e-10	3.9500e-10	4.9400e-10	5.8800e-10
4.1300e-15	3.7800e-12	3.3400e-11	9.7800e-11	1.8500e-10	2.8300e-10	3.3300e-10	3.8300e-10	4.8000e-10	5.7100e-10
2.8100e-15	3.3900e-12	3.1400e-11	9.3500e-11	1.7800e-10	2.7400e-10	3.2200e-10	3.7100e-10	4.6500e-10	5.5500e-10
1.8700e-15	3.0300e-12	2.9500e-11	8.9300e-11	1.7200e-10	2.6400e-10	3.1200e-10	3.5900e-10	4.5100e-10	5.3900e-10
1.2200e-15	2.7100e-12	2.7700e-11	8.5300e-11	1.6500e-10	2.5500e-10	3.0200e-10	3.4800e-10	4.3800e-10	5.2300e-10
7.7300e-16	2.4100e-12	2.6000e-11	8.1400e-11	1.5900e-10	2.4700e-10	2.9200e-10	3.3700e-10	4.2500e-10	5.0800e-10
4.7900e-16	2.1400e-12	2.4400e-11	7.7700e-11	1.5300e-10	2.3800e-10	2.8200e-10	3.2600e-10	4.1200e-10	4.9400e-10
2.8800e-16	1.9000e-12	2.2900e-11	7.4200e-11	1.4700e-10	2.3000e-10	2.7300e-10	3.1600e-10	4.0000e-10	4.7900e-10
1.6800e-16	1.6800e-12	2.1400e-11	7.0800e-11	1.4100e-10	2.2200e-10	2.6400e-10	3.0600e-10	3.8800e-10	4.6500e-10
9.4600e-17	1.4800e-12	2.0100e-11	6.7500e-11	1.3600e-10	2.1500e-10	2.5600e-10	2.9600e-10	3.7600e-10	4.5200e-10
5.1200e-17	1.3000e-12	1.8800e-11	6.4300e-11	1.3100e-10	2.0800e-10	2.4700e-10	2.8700e-10	3.6500e-10	4.3900e-10
2.6600e-17	1.1400e-12	1.7500e-11	6.1300e-11	1.2600e-10	2.0000e-10	2.3900e-10	2.7800e-10	3.5400e-10	4.2600e-10
1.3200e-17	9.9300e-13	1.6400e-11	5.8400e-11	1.2100e-10	1.9400e-10	2.3100e-10	2.6900e-10	3.4300e-10	4.1400e-10
6.1900e-18	8.6400e-13	1.5300e-11	5.5700e-11	1.1600e-10	1.8700e-10	2.2400e-10	2.6000e-10	3.3300e-10	4.0200e-10
2.7400e-18	7.4900e-13	1.4200e-11	5.3000e-11	1.1200e-10	1.8000e-10	2.1600e-10	2.5200e-10	3.2300e-10	3.9000e-10
1.1300e-18	6.4700e-13	1.3200e-11	5.0400e-11	1.0700e-10	1.7400e-10	2.0900e-10	2.4400e-10	3.1300e-10	3.7900e-10
4.3500e-19	5.5700e-13	1.2300e-11	4.8000e-11	1.0300e-10	1.6800e-10	2.0200e-10	2.3600e-10	3.0400e-10	3.6800e-10
1.5300e-19	4.7700e-13	1.1500e-11	4.5700e-11	9.9000e-11	1.6200e-10	1.9500e-10	2.2900e-10	2.9400e-10	3.5700e-10
1.3100e-22	2.0600e-13	7.8300e-12	3.5400e-11	8.0700e-11	1.3600e-10	1.6500e-10	1.9400e-10	2.5200e-10	3.0800e-10
1.9600e-42	2.4300e-14	3.3600e-12	2.0600e-11	5.2900e-11	9.4300e-11	1.1700e-10	1.3900e-10	1.8500e-10	2.2800e-10
9.2300e-145	5.9300e-15	2.0800e-12	1.5400e-11	4.2500e-11	7.8300e-11	9.7800e-11	1.1800e-10	1.5800e-10	1.9600e-10
1.1915e-06	1.0200e-15	1.2300e-12	1.1400e-11	3.3900e-11	6.4800e-11	8.1800e-11	9.9300e-11	1.3500e-10	1.6900e-10];


    De_est = interp2(T_sample,c_e_sample,De_sample,T,c_e,'linear');
end
