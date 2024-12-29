param = getparam('P4-2v');

% Add transmit apodization

param.TXapodization = cos(linspace(-1,1,64)*pi/2);

% set delays
tilt = deg2rad(linspace(-20,20,7)); % tilt angles in rad
txdel = cell(7,1); % this cell will contain the transmit delays
for k = 1:7
    txdel{k} = txdelay(param,tilt(k),deg2rad(60));
end
% Plot the transmit delays
stem(txdel{1}*1e6)  % Transmit delay in microseconds
xlabel('Element number')  % Label for x-axis
ylabel('Delays (\mus)')  % Label for y-axis (microseconds)
title('TX delays for a 60{\circ}-wide -20{\circ}-tilted wave')  % Title of the plot
axis tight square  % Adjust the axis and set aspect ratio to be square

% Save the figure
saveas(gcf, 'TX_delays_plot.png');  % Save the plot as a PNG file

%Define a 100 $\times$ 100 polar grid using IMPOLGRID:

[xi,zi] = impolgrid([100 100],15e-2,deg2rad(120),param);
option.WaitBar = true;
P = pfield(xi,0*xi,zi,txdel{1},param,option);
% Plot the RMS pressure field
pcolor(xi*1e2, zi*1e2, 20*log10(P/max(P(:))))  % Pressure field in dB
shading interp  % Interpolated shading
xlabel('x (cm)')  % Label for x-axis
ylabel('z (cm)')  % Label for z-axis
title('RMS pressure field for a 60{\circ}-wide -20{\circ}-tilted wave')  % Title of the plot
axis equal ij tight  % Adjust axis and aspect ratio
caxis([-20 0])  % Set the dynamic range of the color axis
cb = colorbar;  % Add a colorbar
cb.YTickLabel{end} = '0 dB';  % Set the last colorbar label to '0 dB'
colormap(hot)  % Use the hot colormap

% Save the figure
saveas(gcf, 'RMS_pressure_field.png');  % Save the plot as a PNG file

I = rgb2gray(imread('heart.jpg'));
% Pseudorandom distribution of scatterers (depth is 15 cm)
[x,y,z,RC] = genscat([NaN 15e-2],1540/param.fc,I);
% Plot the scatterers for the cardiac 5-chamber view
scatter(x*1e2, z*1e2, 2, abs(RC).^.25, 'filled')  % Scatter plot with color based on reflection coefficients
colormap([1-hot; hot])  % Use a colormap that transitions from white to hot
axis equal ij tight  % Set equal scaling and reverse the z-axis
set(gca, 'XColor', 'none', 'box', 'off')  % Remove x-axis color and box
title('Scatterers for a cardiac 5-chamber view')  % Set the plot title
ylabel('[cm]')  % Label the y-axis with units in centimeters

% Save the figure
saveas(gcf, 'Cardiac_scatterers.png');  % Save the scatter plot as a PNG file

% Initialize a cell array to store the RF signals for each of the 7 series
RF = cell(7,1); 

% Set the sampling frequency for RF signals, based on the center frequency
param.fs = 4 * param.fc;  % Sampling frequency in Hz (4 times the center frequency)

% Remove the waitbar during simulation for faster execution
option.WaitBar = false; 

% Initialize a waitbar to display progress of the RF signal simulation
h = waitbar(0, ''); 

% Loop through the 7 tilt angles to simulate RF signals for each tilt angle
for k = 1:7
    % Update the waitbar to show progress for the current series
    waitbar(k/7, h, ['SIMUS: RF series #' int2str(k) ' of 7'])
    
    % Call the SIMUS function to simulate RF signals based on the scatterer data, transmit delays, and other parameters
    RF{k} = simus(x, y, z, RC, txdel{k}, param, option);
end

% Close the waitbar once the simulation is complete
close(h);

% Extract the 32nd RF signal from the first series (index 1)
rf = RF{1}(:, 32);

% Create a time vector (in microseconds), assuming the sampling frequency is 'param.fs' in Hz
t = (0:numel(rf)-1) / param.fs * 1e6;  % time in microseconds

% Plot the RF signal against time
figure;
plot(t, rf);

% Customize the axes and appearance
set(gca, 'YColor', 'none', 'box', 'off');  % Remove the Y-axis ticks and box border
xlabel('Time (\mus)');  % Label for X-axis (time in microseconds)
title('RF signal of the 32^{nd} element (1^{st} series, tilt = -20{\circ})');  % Title of the plot
axis tight;  % Adjust axis limits to fit the data tightly

saveas(gcf, 'RF_signal_32_element.png');

IQ = cell(7,1);  % this cell will contain the I/Q series

for k = 1:7
    IQ{k} = rf2iq(RF{k},param.fs,param.fc);
end
iq = IQ{1}(:,32);
plot(t,real(iq),t,imag(iq))
set(gca,'YColor','none','box','off')
xlabel('time (\mus)')
title('I/Q signal of the 32^{th} element (1^{st} series, tilt = -20{\circ})')
legend({'in-phase','quadrature'})
axis tight
saveas(gcf, 'IQ_signal_32_element.png');

[xi,zi] = impolgrid([256 128],15e-2,deg2rad(80),param);

bIQ = zeros(256,128,7);  % this array will contain the 7 I/Q images

h = waitbar(0,'');
for k = 1:7
    waitbar(k/7,h,['DAS: I/Q series #' int2str(k) ' of 7'])
    bIQ(:,:,k) = das(IQ{k},xi,zi,txdel{k},param);
end
close(h)
bIQ = tgc(bIQ);

I = bmode(bIQ(:,:,1),50); % log-compressed image
pcolor(xi*1e2,zi*1e2,I)
shading interp, colormap gray
title('DW-based echo image with a tilt angle of -20{\circ}')

axis equal ij
set(gca,'XColor','none','box','off')
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-50 dB','0 dB'};
ylabel('[cm]')
cIQ = sum(bIQ,3); % this is the compound beamformed I/Q
I = bmode(cIQ,50); % log-compressed image
pcolor(xi*1e2,zi*1e2,I)
shading interp, colormap gray
title('Compound DW-based cardiac echo image')

axis equal ij
set(gca,'XColor','none','box','off')
c = colorbar;
c.YTick = [0 255];
c.YTickLabel = {'-50 dB','0 dB'};
ylabel('[cm]')
