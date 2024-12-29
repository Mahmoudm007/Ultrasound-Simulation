% Load phased-array transducer parameters
param = getparam('P4-2v'); % Example transducer: 64-element cardiac phased array

% Define common simulation parameters
param.fs = 4 * param.fc;  % Sampling frequency
param.c = 1540;           % Speed of sound [m/s]
depth = 10e-2;            % Depth of simulation area [m]
width = 10e-2;            % Width of simulation area [m]
grid_points = [256, 256]; % Number of grid points for pressure field visualization

%% 1. Focused Pressure Field (Single Focus)
focus_depth = 6e-2;       % Focal point [m]
txdel_focused = txdelay(param, focus_depth); % Calculate transmit delays
[x_focused, z_focused] = impolgrid(grid_points, depth, pi/3, param);
p_focused = pfield(x_focused, z_focused, txdel_focused, param);

% Display focused field
figure;
pcolor(x_focused * 100, z_focused * 100, p_focused);
shading interp; colormap hot; colorbar;
title('Focused Pressure Field');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;

%% 2. Diverging Pressure Field
txdel_diverging = txdelay(param, 0, pi/3); % Generate diverging wave delays
[x_diverging, z_diverging] = impolgrid(grid_points, depth, pi/3, param);
p_diverging = pfield(x_diverging, z_diverging, txdel_diverging, param);

% Display diverging field
figure;
pcolor(x_diverging * 100, z_diverging * 100, p_diverging);
shading interp; colormap hot; colorbar;
title('Diverging Pressure Field');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;

%% 3. Multi-Focus Pressure Field
focus_depths = [4e-2, 6e-2, 8e-2]; % Multiple focal points [m]
p_multifocus = zeros(size(x_focused));
for focus = focus_depths
    txdel_multi = txdelay(param, focus); % Transmit delays for each focal point
    p_multifocus = p_multifocus + pfield(x_focused, z_focused, txdel_multi, param);
end

% Display multi-focus field
figure;
pcolor(x_focused * 100, z_focused * 100, p_multifocus);
shading interp; colormap hot; colorbar;
title('Multi-Focus Pressure Field');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;

%% 4. Focused Pressure Field on a Polar Grid
% Define a focus position
xf = 2e-2; 
zf = 5e-2;

% Calculate the transmit delays with TXDELAY
txdel = txdelay(xf, zf, param);

% Create a 60-degree wide, 10-cm deep, 100x50 polar grid
[x_polar, z_polar] = impolgrid([100 50], 10e-2, pi/3, param);

% Simulate the RMS pressure field with PFIELD
y_polar = zeros(size(x_polar)); 
P_polar = pfield(x_polar, y_polar, z_polar, txdel, param);

% Display a scatter plot of the pressure field
figure;
scatter(x_polar(:) * 1e2, z_polar(:) * 1e2, 5, 20 * log10(P_polar(:) / max(P_polar(:))), 'filled');
colormap jet; 
caxis([-20 0]);
colorbar; 
axis equal ij tight;
xlabel('[cm]'); 
ylabel('[cm]');
title('Scatter Plot of Pressure Field');

% Display the pressure field
figure;
pcolor(x_polar * 1e2, z_polar * 1e2, 20 * log10(P_polar / max(P_polar(:))));
shading interp; 
colormap hot; 
axis equal ij tight;
xlabel('[cm]');
ylabel('[cm]');
caxis([-20 0]);
colorbar;
title('Pressure Field on Polar Grid');

%% 5. Diverging Wave Backpropagation
% Design scatterers for the phased-array simulation
xs = [4.3 3.2 1.7 0 -1.7 -3.2 -4.3 0 -2.5 2.5] * 1e-2; % Scatterer x-positions [m]
zs = [7 8.1 8.8 9 8.8 8.1 7 5 2 2] * 1e-2;             % Scatterer z-positions [m]

% Set all transmit delays to 0
txdel_backprop = zeros(1, param.Nelements);

% Calculate and display the pressure field
[x_backprop, z_backprop] = meshgrid(linspace(-.1, .1, 256), linspace(0, .1, 128));
P_backprop = pfield(x_backprop, z_backprop, txdel_backprop, param);

% Display the backpropagated pressure field
figure;
imagesc(x_backprop(1,:) * 1e2, z_backprop(:,1) * 1e2, 20 * log10(P_backprop / max(P_backprop, [], 'all')));
axis equal ij tight;
caxis([-20 0]); % Dynamic range = [-20,0] dB
colorbar;
colormap hot;
xlabel('x (cm)');
ylabel('z (cm)');
title('Diverging Wave from Phased Array');
hold on;
plot(xs * 1e2, zs * 1e2, 'bo', 'LineWidth', 2);
legend('Scatterers');
hold off;

% Add a "ghost" scatterer to ensure RF is simulated to 10-cm depth
RC = [ones(1, numel(xs)), 0]; 
xs = [xs, 0]; 
zs = [zs, 1e-2];

% Create a GIF movie named 'phasedArraySmiley.gif'
param.movie = [10 10 30]; % 10-cm-by-10-cm ROI, resolution = 30 pix/cm
mkmovie(xs, zs, RC, txdel_backprop, param, 'phasedArraySmiley.gif');

%% 6. 2-D Focused Pressure Field with a Phased-Array Transducer
xf = 2e-2; zf = 5e-2; % Focus position [m]
txdel = txdelay(xf, zf, param);
x = linspace(-4e-2, 4e-2, 200); % X grid [m]
z = linspace(0, 10e-2, 200);    % Z grid [m]
[x, z] = meshgrid(x, z);
P = pfield(x, z, txdel, param);

% Display the pressure field
figure;
imagesc(x(1,:) * 1e2, z(:,1) * 1e2, 20 * log10(P / max(P(:))));
caxis([-20 0]);
colorbar;
colormap hot;
axis equal ij tight;
xlabel('x (cm)');
ylabel('z (cm)');
title('2-D Focused Pressure Field');
hold on;
plot(xf * 1e2, zf * 1e2, 'bo', 'MarkerFaceColor', 'b');
legend('Focus Point', 'Location', 'South');
hold off;

%% 7. 3-D Focused Pressure Field with a Phased-Array Transducer
xf = -2e-2; zf = 5e-2; % Focus position [m]
txdel = txdelay(xf, zf, param);
x = linspace(-4e-2, 4e-2, 200); % X grid [m]
y = linspace(-0.75, 0.75, 50) * param.height; % Y grid [m]
z = linspace(0, 10e-2, 200);  % Z grid [m]

% Azimuthal plane
[xaz, zaz] = meshgrid(x, z);
yaz = zeros(size(xaz));
Paz = pfield(xaz, yaz, zaz, txdel, param);

% Elevation plane
[yel, zel] = meshgrid(y, z);
xel = ones(size(yel)) * xf;
Pel = pfield(xel, yel, zel, txdel, param);

% Focal plane
[xfo, yfo] = meshgrid(x, y);
zfo = ones(size(xfo)) * zf;
Pfo = pfield(xfo, yfo, zfo, txdel, param);

% Display the 3-D pressure field
figure;
Pmax = max(Paz(:));
surf(xaz, yaz, zaz, 20 * log10(Paz / Pmax), 'EdgeColor', 'none');
hold on;
surf(xel, yel, zel, 20 * log10(Pel / Pmax), 'EdgeColor', 'none');
surf(xfo, yfo, zfo, 20 * log10(Pfo / Pmax), 'EdgeColor', 'none');
hold off;
axis equal;
set(gca, 'ZDir', 'reverse');
title('3-D Focused Pressure Field');
caxis([-20 0]);
colormap hot;

%% Diverging wave with a phased array
param = getparam('P4-2v');

width = 75/180*pi; % width angle in rad 
tilt = 0; % tilt angle in rad 
txdel = txdelay(param,tilt,width); % in s

param.TXapodization = cos(linspace(-3*pi/8,3*pi/8,64));
bar(param.TXapodization)

xlabel('Element number') 
title('TX apodization') 
axis tight square

x = linspace(-4e-2,4e-2,200); % in m
z = linspace(0,10e-2,200); % in m
[x,z] = meshgrid(x,z);

P = pfield(x,z,txdel,param);

imagesc(x(1,:)*1e2,z(:,1)*1e2,20*log10(P/max(P,[],'all')))
caxis([-20 0]) % dynamic range = [-20,0] dB
c = colorbar;
c.YTickLabel{end} = '0 dB';
colormap hot
axis equal ij tight
xlabel('x (cm)'), ylabel('z (cm)')
title('Diverging wave - RMS pressure field')

%% Multi-line transmit with phased-array

x0 = 2e-2; z0 = 5e-2; % in m
xf = [-x0 0 x0]; zf = [z0 sqrt(x0^2+z0^2) z0]; % focus points (in m)
txdel = txdelay(xf,zf,param); % in s

% Define the image grid.
x = linspace(-5e-2,5e-2,200); % in m
z = linspace(0,10e-2,200); % in m
[x,z] = meshgrid(x,z); % image grid
y = zeros(size(x));

P = pfield(x,y,z,txdel,param);

% Display P field
imagesc(x(1,:)*1e2,z(:,1)*1e2,20*log10(P/max(P,[],'all')))
caxis([-20 0]) % dynamic range = [-20,0] dB
c = colorbar;
c.YTickLabel{end} = '0 dB';
colormap hot
axis equal ij tight
xlabel('x (cm)'), ylabel('z (cm)')
title('MLT - RMS pressure field')
hold on

plot(xf*1e2,zf*1e2,'bo','MarkerFaceColor','b')
legend('focus points','Location','South') 
hold off