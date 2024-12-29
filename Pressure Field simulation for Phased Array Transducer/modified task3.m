% Load phased-array transducer parameters
param = getparam('P4-2v'); % Example transducer: 64-element cardiac phased array

% Define common simulation parameters
param.fs = 4 * param.fc;  % Sampling frequency
param.c = 1540;           % Speed of sound [m/s]
depth = 10e-2;            % Depth of simulation area [m]
width = 10e-2;            % Width of simulation area [m]
grid_points = [256, 256]; % Number of grid points for pressure field visualization











1. Focused Pressure Field (Single Focus)
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












6. 2-D Focused Pressure Field with a Phased-Array Transducer
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


































2. Diverging Pressure Field
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














3. Multi-Focus Pressure Field
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