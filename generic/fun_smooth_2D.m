function [smooth_P, xp_new, yp_new] = fun_smooth_2D(xp,yp,P,C)

% Define your matrices C and P here
% C = ...
% P = ...



% Obtain the dimensions of the matrices
[m, n] = size(P);

% Create a weighted version of P
weightedP = P .* C;

% Generate the original grid
[x, y] = meshgrid(1:n, 1:m);

% Flatten the original grid and weighted P

% Flatten the original grid and weighted P
x_flat = x(:);
y_flat = y(:);
weightedP_flat = weightedP(:);

% Define the new resolution you want for the output grid
resolution_multiplier = 5; % Increase this value for higher resolution

xp_new = linspace(xp(1), xp(end), n * resolution_multiplier);
yp_new = linspace(yp(1), yp(end), n * resolution_multiplier);

% Generate the new, higher resolution grid
[x_new, y_new] = meshgrid(linspace(1, n, n * resolution_multiplier), linspace(1, m, m * resolution_multiplier));

% Interpolate the weighted P values onto the new grid using griddata
smoothed_weightedP = griddata(x_flat, y_flat, weightedP_flat, x_new, y_new, 'cubic');

% Normalize the smoothed weighted P by interpolating the weights (C matrix) onto the new grid
C_flat = C(:);
smoothed_C = griddata(x_flat, y_flat, C_flat, x_new, y_new, 'cubic');

% Take care of any NaN values that might appear due to interpolation
smoothed_C(isnan(smoothed_C)) = 0;
smoothed_weightedP(isnan(smoothed_weightedP)) = 0;

% Compute the smooth P function
smooth_P = smoothed_weightedP ./ smoothed_C;
smooth_P(isnan(smooth_P)) = 0; % Set any NaN values to 0

% Plot the original P matrix and the smooth P function for comparison
figure;
subplot(1, 2, 1);
imagesc(P);
title('Original P Matrix');
xlabel('n');
ylabel('m');
colorbar;

subplot(1, 2, 2);
imagesc(smooth_P);
title('Smooth P Function');
xlabel('n (Higher Resolution)');
ylabel('m (Higher Resolution)');
colorbar;