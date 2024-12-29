% Load linear-array transducer parameters
param = getparam('L11-5v'); % Example transducer: 128-element linear array transducer

% Define common simulation parameters
param.fs = 4 * param.fc;  % Sampling frequency
param.c = 1540;           % Speed of sound [m/s]
depth = 10e-2;            % Depth of simulation area [m]
width = 10e-2;            % Width of simulation area [m]
grid_points = [256, 256]; % Number of grid points for pressure field visualization

%% 1. Focused Pressure Field (Single Focus)
focus_depth = 6e-2;       % Focal point [m]
txdel_focused = txdelay(param, focus_depth); % Calculate transmit delays
[x_focused, z_focused] = meshgrid(linspace(-width/2, width/2, grid_points(1)), ...
                                  linspace(0, depth, grid_points(2)));
p_focused = pfield(x_focused, z_focused, txdel_focused, param);

% Display focused field
figure;
pcolor(x_focused * 100, z_focused * 100, 20 * log10(p_focused / max(p_focused, [], 'all')));
shading interp; colormap hot; colorbar;
caxis([-20 0]); % Dynamic range: [-20, 0] dB
title('Focused Pressure Field');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;

%% 2. Plane Wave
tilt_angle = 0; % No tilt for a vertical plane wave
txdel_plane = txdelay(param, tilt_angle); % Generate plane wave delays
[x_plane, z_plane] = meshgrid(linspace(-width/2, width/2, grid_points(1)), ...
                              linspace(0, depth, grid_points(2)));
p_plane = pfield(x_plane, z_plane, txdel_plane, param);

% Display plane wave field
figure;
pcolor(x_plane * 100, z_plane * 100, 20 * log10(p_plane / max(p_plane, [], 'all')));
shading interp; colormap hot; colorbar;
caxis([-20 0]); % Dynamic range: [-20, 0] dB
title('Plane Wave Pressure Field');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;

%% 3. Multi-Focus Pressure Field
focus_depths = [4e-2, 6e-2, 8e-2]; % Multiple focal depths [m]
p_multifocus = zeros(size(x_focused));
for focus = focus_depths
    txdel_multi = txdelay(param, focus); % Transmit delays for each focus
    p_multifocus = p_multifocus + pfield(x_focused, z_focused, txdel_multi, param);
end

% Display multi-focus field
figure;
pcolor(x_focused * 100, z_focused * 100, 20 * log10(p_multifocus / max(p_multifocus, [], 'all')));
shading interp; colormap hot; colorbar;
caxis([-20 0]); % Dynamic range: [-20, 0] dB
title('Multi-Focus Pressure Field');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;

%% 4. Plane Wave with Apodization
apod = hanning(param.Nelements)'; % Example: Hanning window apodization
param.TXapodization = apod; % Apply the apodization to the transducer
txdel_plane_apod = txdelay(param, tilt_angle); % Plane wave with apodization

% Simulate pressure field with apodization
p_plane_apod = pfield(x_plane, z_plane, txdel_plane_apod, param);

% Display plane wave with apodization
figure;
pcolor(x_plane * 100, z_plane * 100, 20 * log10(p_plane_apod / max(p_plane_apod, [], 'all')));
shading interp; colormap hot; colorbar;
caxis([-20 0]); % Dynamic range: [-20, 0] dB
title('Plane Wave with Apodization');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;

%% 5. Multi-Line Transmit (MLT)
focus_offsets = [-2e-2, 0, 2e-2]; % Three focal positions (lateral)
focus_depths = [5e-2, 5e-2, 5e-2]; % All at the same depth
txdel_mlt = txdelay(focus_offsets, focus_depths, param); % Calculate MLT delays

% Simulate MLT pressure field
p_mlt = pfield(x_focused, z_focused, txdel_mlt, param);

% Display MLT pressure field
figure;
pcolor(x_focused * 100, z_focused * 100, 20 * log10(p_mlt / max(p_mlt, [], 'all')));
shading interp; colormap hot; colorbar;
caxis([-20 0]); % Dynamic range: [-20, 0] dB
title('Multi-Line Transmit (MLT) Pressure Field');
xlabel('Lateral Position [cm]');
ylabel('Depth [cm]');
axis equal tight;
hold on;
plot(focus_offsets * 100, focus_depths * 100, 'bo', 'MarkerFaceColor', 'b');
legend('Focus Points', 'Location', 'South');
hold off;

%% Focused pressure field with a subaperture of a linear array
L = param.pitch*(param.Nelements-1); % array aperture (in m)
xe = linspace(-0.5,0.5,param.Nelements)*L;
ze = zeros(1,param.Nelements);

param_suba = param;
param_suba.Nelements = 24;

xf = 0; zf = 2.5e-2; % focus position (in m)

txdel_suba = txdelay(xf,zf,param_suba); % in s

txdel_1to24 = NaN(1,128); 
txdel_1to24(1:24) = txdel_suba;

txdel_91to114 = NaN(1,128);
txdel_91to114(91:114) = txdel_suba;

x = linspace(-2.5e-2,2.5e-2,150); % in m
z = linspace(0,5e-2,150); % in m
[x,z] = meshgrid(x,z);

P_1to24 = pfield(x,z,txdel_1to24,param);
P_91to114 = pfield(x,z,txdel_91to114,param);

% Display the acoustic pressure field when firing elements #1 to #24.
imagesc(x(1,:)*1e2,z(:,1)*1e2,20*log10(P_1to24/max(P_1to24,[],'all')))
caxis([-20 0]) % dynamic range = [-20,0] dB
c = colorbar;
c.YTickLabel{end} = '0 dB';
colormap hot
axis equal ij tight
xlabel('x (cm)'), ylabel('z (cm)')
title('Focused wave with elements #1 to #24')

hold on
plot(xe*1e2,ze*1e2,'g','Linewidth',5)
plot(xe*1e2,ze*1e2,'g','Linewidth',5)
hold off

%% Display the acoustic pressure field when firing elements #91 to #114.
imagesc(x(1,:)*1e2,z(:,1)*1e2,20*log10(P_91to114/max(P_91to114,[],'all')))

caxis([-20 0]) % dynamic range = [-20,0] dB
c = colorbar;
c.YTickLabel{end} = '0 dB';
colormap hot
axis equal ij tight
xlabel('x (cm)'), ylabel('z (cm)')
title('Focused wave with elements #91 to #114')

hold on
plot(xe*1e2,ze*1e2,'g','Linewidth',5)
plot(xe(91:114)*1e2,ze(91:114)*1e2,'r','Linewidth',5)
hold off

%% Steered plane wave with a linear array
tilt = 10/180*pi; % tilt angle in rad 
txdel = txdelay(param,tilt); % in s

x = linspace(-4e-2,4e-2,150); % in m 
z = linspace(0,8e-2,150); % in m 
[x,z] = meshgrid(x,z);

P = pfield(x,z,txdel,param);


imagesc(x(1,:)*1e2,z(:,1)*1e2,20*log10(P/max(P,[],'all'))) 
caxis([-20 0]) % dynamic range = [-20,0] dB 
c = colorbar; 
c.YTickLabel{end} = '0 dB'; 
colormap hot 

axis equal ij tight 
xlabel('x (cm)'), ylabel('z (cm)') 
title('A steered plane wave with a linear array') 
% Calculate the positions of the element centers. 
L = param.pitch*(param.Nelements-1); % array aperture (in m) 
xe = linspace(-0.5,0.5,param.Nelements)*L; 
ze = zeros(1,param.Nelements); 
hold on 
plot(xe*1e2,ze*1e2,'g','Linewidth',5) 
hold off

%%
