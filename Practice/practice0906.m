

% frequency vector
f = logspace(-3, 3, 120);
w = 2*pi*f;

% parameter
A=0.1;
R=1;
C=1;
params = [R, C, A];

% model calculation
z_mod = z_model(w,params);


% noise
z_re_noise = (3*R/100)*(rand(size(z_mod))-0.5);
z_im_noise = (4*R/100)*(rand(size(z_mod))-0.5);
z_syn = z_mod +z_re_noise + z_im_noise*1i;
% improve: use normal distribution randome number with mean and std

% plot
figure(1)
plot(real(z_mod), -imag(z_mod),'linewidth',2)
hold on
plot(real(z_syn),-imag(z_syn),'o','markersize',4,'linewidth',0.5)

% fitting
% discuss with the team


% plot by the fitted parameter

function [cost] = rmse(z_model, z_data)

cost = sqrt(sum((real(z_model - z_data)).^2 +(imag(z_model - z_data)).^2));


end



function [Z] = z_model(w,params)
R=params(1);
C=params(2);
A=params(3);

Z_W = A .* (1 - 1i) ./ sqrt(w);
Z_RW = R + Z_W;
Z_C = 1 ./ (1i*w*C);
Z = (Z_RW .* Z_C) ./ (Z_RW + Z_C);
end

