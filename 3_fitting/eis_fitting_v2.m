
% This function (1) finds the best-fit parameters for a given full cell eis spectra. 
%               (2) plot the optimized solution and the experimenttal data.
% Calls a function "godfrey_JS_V6" that calculates the full cell impedance.
% [v2]  taking V6b function, which has other outputs
%       giving confidence interval.


%% 0. Condition at which the model is optimized
    soc0 = 48; % [%]
    T0 = 298.15; % [K]

%% 1. Read the Experimental Specta to Fit
    exp_data = importdata('C:\Users\jsong\Documents\MATLAB\BSL_EIS\4_data\example1_eis_fc_soc48_t25.txt');
    w_exp = exp_data.data(:,1);         % [Hz]
    z_exp_real = exp_data.data(:,2);    % [Ohm]
    z_exp_imag = - exp_data.data(:,3);    % [Ohm]
    z_exp_realmat = [z_exp_real z_exp_imag];    % [Ohm]

%% 2. Define the Weighting vector
    w_switch = 0; % choose between minimizing (0) absolute error; (1) relative error
    % if minimizing the relative error
    if w_switch == 1
    weight = (z_exp_real.^2 + z_exp_imag.^2).^0.5;
    weight_matrix = [weight weight];
    % if minimizing the absolute error
    elseif w_switch == 0
    weight_matrix = ones(size(z_exp_realmat));
    end
%% 3. Initial Values of the Fit Parameters
    factor0(1) = 1; % Factor for external circuit resistance
    factor0(2) = 1; % Factor for Bruggman alpha
    factor0(3) = 2; % Factor for p_am1_Ds
    factor0(4) = 10; % Factor for n_am1_Ds
    factor0(5) = 1; % Factor for p_am1_i0
    factor0(6) = 10; % Factor for p_Cdl
    factor0(7) = 1; % Factor for p_am1_dUdc
    factor0(8) = 1; % Factor for n_am1_dUdc
    factor0(9) = 1; % Factor for n_Cdl
    factor0(10) = 0.5; % Factor for n_am1_i0
    
%% 4. Perform the Optimization
    % optimization options
    options = optimset('display','iter','MaxIter',400,'MaxFunEvals',1e5,...
        'TolFun',1e-10,'TolX',1e-8,'FinDiffType','central');
    % bounds: [ low   upper ]           
    bounds = [  -1,   1.5;      % R_itsc 
                0.65,  1.6;      % burg
                0.1,    20;     % p_Ds
                0.1,    20;     % n_Ds
                0.5,    1.5;    % p_i0
                0.01,   1000;   % p_Cdl
                0.8,    1.2;    % p_dUdc
                0.1,    10;    % n_dUdc
                0.01,   1000;  % n_Cdl
                0.5,    20];     %n_i0
    lb = bounds(:,1).'; ub = bounds(:,2).';
    % call the impedance function
    eis_model_weighted = @(factors,w)JS_EIS_model_V6b(w,factors,soc0,T0).*weight_matrix;
    z_exp_realmat_weighted = z_exp_realmat.*weight_matrix; 
    tic;
    [factor_hat,resnorm,residual,~,~,~,jacobian_hat] = lsqcurvefit(eis_model_weighted,factor0,w_exp,z_exp_realmat_weighted,lb,ub,options);
    toc;
    conf_vec = nlparci(factor_hat,residual,'Jacobian',jacobian_hat);
%% 5. Plot and Compare
    [z_mod0_realmat,used_para0,z_local0] = JS_EIS_model_V6b(w_exp,factor0,soc0,T0);
    [z_mod1_realmat,used_para1,z_local1] = JS_EIS_model_V6b(w_exp,factor_hat,soc0,T0);
    
    
    % Nyquist plot
    figure(1); hold on;
    plot(z_exp_realmat(:,1),-z_exp_realmat(:,2),'ok','linewidth',1)
    plot(z_mod0_realmat(:,1),-z_mod0_realmat(:,2),'ob','linewidth',1)
    plot(z_mod1_realmat(:,1),-z_mod1_realmat(:,2),'or','linewidth',1)
    legend('experimental','before-fitting','after-fitting')
    axis_limit = 1.1*max(max(abs(z_exp_realmat)));
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman',... % Fonts
    'XLim',[0 axis_limit],'Ylim',[0 axis_limit])
    hold off
    
    % Nyquist plot for local impedance
    figure(2); hold on;
    plot(real(z_local0(:,2)),-imag(z_local0(:,2)),'--g',real(z_local1(:,2)),-imag(z_local1(:,2)),'-g')
    plot(real(z_local0(:,3)),-imag(z_local0(:,3)),'--r',real(z_local1(:,3)),-imag(z_local1(:,3)),'-r')
    plot(real(z_local0(:,4)),-imag(z_local0(:,4)),'--b',real(z_local1(:,4)),-imag(z_local1(:,4)),'-b')
    plot(real(z_local0(:,5)),-imag(z_local0(:,5)),'--c',real(z_local1(:,5)),-imag(z_local1(:,5)),'-c')
    legend('fullcell0','fullcell1','cathode0','cathode1','anode0','anode1','separator0','separator1')
    set(gca,'Box','on',... %Axis Properties: BOX   
    'PlotBoxAspectRatio',[1 1 1],... % Size - you can either use 'position' or 'dataaspectratio' or their combinations
    'FontUnits','points','FontSize',10,'FontName','Times New Roman')
    hold off
    
    % Bode plots
    % skipped
