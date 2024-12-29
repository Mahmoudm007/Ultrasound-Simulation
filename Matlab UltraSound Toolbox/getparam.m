function param = getparam(probe)

%GETPARAM   Get parameters of a uniform linear or convex array
%   PARAM = GETPARAM opens a dialog box which allows you to select a
%   transducer whose parameters are returned in PARAM.
%
%   PARAM = GETPARAM(PROBE), where PROBE is a string, returns the prameters
%   of the transducer given by PROBE.
%
%   The structure PARAM is used in several functions of MUST (Matlab
%   UltraSound Toolbox). The structure returned by GETPARAM contains only
%   the fields that describe a transducer. Other fields may be required in
%   some MUST functions.
%
%   PROBE can be one of the following:
%   ---------------------------------
%     1) 'L11-5v' (128-element, 7.6-MHz linear array)
%     2) 'L12-3v' (192-element, 7.5-MHz linear array)
%     3) 'C5-2v' (128-element, 3.6-MHz convex array)
%     4) 'P4-2v' (64-element, 2.7-MHz phased array)
%     5) 'Vermon-3' (1024-element, 3-MHz matrix array)
%     6) 'Vermon-8' (1024-element, 8-MHz matrix array)
%
%   These are <a
%   href="matlab:web('https://verasonics.com/verasonics-transducers/')">Verasonics</a> and <a href="matlab:web('https://verasonics.com/matrix-array/')">Vermon</a> transducers.
%   Feel free to complete this list for your own use.
%
%   PARAM is a structure that contains the following fields:
%   --------------------------------------------------------
%      For PFIELD and SIMUS:
%      --------------------
%   1) PARAM.Nelements: number of elements in the transducer array
%   2) PARAM.fc: center frequency (in Hz)
%   3) PARAM.pitch: element pitch (in m)
%   4) PARAM.width: element width (in m)
%   5) PARAM.kerf: kerf width (in m)
%   6) PARAM.bandwidth: 6-dB fractional bandwidth (in %)
%   7) PARAM.radius: radius of curvature (in m, Inf for a linear array)
%   8) PARAM.focus: elevation focus (in m)
%   9) PARAM.height: element height (in m)
%
%      For PFIELD3 and SIMUS3:
%      ----------------------
%   1) PARAM.Nelements: number of elements in the transducer array
%   2) PARAM.fc: center frequency (in Hz)
%   3) PARAM.width: element width, in the x-direction (in m)
%   4) PARAM.height: element height, in the y-direction (in m)
%   5) PARAM.bandwidth: 6-dB fractional bandwidth (in %)
%   7) PARAM.radius: radius of curvature (in m, Inf for a linear array)
%   8) PARAM.elements: x- and y-coordinates of the element centers (in m).
%      It is a two-row matrix, with the 1st and 2nd rows containing the x
%      and y coordinates, respectively.
%
%
%   Example:
%   -------
%   %-- Generate a focused pressure field with a phased-array transducer
%   % Phased-array @ 2.7 MHz:
%   param = getparam('P4-2v');
%   % Focus position:
%   x0 = 2e-2; z0 = 5e-2;
%   % TX time delays:
%   dels = txdelay(x0,z0,param);
%   % Grid:
%   x = linspace(-4e-2,4e-2,200);
%   z = linspace(param.pitch,10e-2,200);
%   [x,z] = meshgrid(x,z);
%   y = zeros(size(x));
%   % RMS pressure field:
%   P = pfield(x,y,z,dels,param);
%   imagesc(x(1,:)*1e2,z(:,1)*1e2,20*log10(P/max(P(:))))
%   hold on, plot(x0*1e2,z0*1e2,'k*'), hold off
%   colormap hot, axis equal tight
%   caxis([-20 0])
%   c = colorbar;
%   c.YTickLabel{end} = '0 dB';
%   xlabel('[cm]')
%
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also TXDELAY, TXDELAY3, PFIELD, PFIELD3, SIMUS, SIMUS3, GETPULSE.
%
%   -- Damien Garcia -- 2015/03, last update: 2022/12
%   -- updated by François Varray (Vermon matrix arrays) -- 2022/11
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>

if nargin==0
    ProbeList = {'L11-5v','L12-3v','C5-2v','P4-2v',...
        'Vermon-3','Vermon-8'};
    ProbeList = sort(ProbeList);
    selection = listdlg('ListString',ProbeList,...
        'PromptString','Select a probe:',...
        'SelectionMode','single');
    try
        probe = ProbeList{selection};
    catch
        param = struct([]);
        return
    end
end

% from computeTrans.m (Verasonics, version post Aug 2019)
switch upper(probe)
    case 'L11-5V'
        % --- L11-5v (Verasonics) ---
        param.fc = 7.6e6; % Transducer center frequency [Hz]
        param.kerf = 30e-6; % Kerf [m]
        param.width = 270e-6; % Width [m]
        param.pitch = 300e-6; % Pitch [m]
        param.Nelements = 128;
        param.bandwidth = 77; % Fractional bandwidth [%]
        param.radius = Inf;
        param.height = 5e-3; % Elevation height [m]
        param.focus = 18e-3; % Elevation focus [m]
    case 'L12-3V'
        % --- L12-3v (Verasonics) ---
        param.fc = 7.54e6; % Transducer center frequency [Hz]
        param.kerf = 30e-6; % Kerf [m]
        param.width = 170e-6; % Width [m]
        param.pitch = 200e-6; % Pitch [m]
        param.Nelements = 192;
        param.bandwidth = 93; % Fractional bandwidth [%]
        param.radius = Inf;
        param.height = 5e-3; % Elevation height [m]
        param.focus = 20e-3; % Elevation focus [m]
    case 'C5-2V'
        % --- C5-2v (Verasonics) ---
        param.fc = 3.57e6; % Transducer center frequency [Hz]
        param.kerf = 48e-6; % Kerf [m]
        param.width = 460e-6; % Width [m]
        param.pitch = 508e-6; % Pitch [m]
        param.Nelements = 128;
        param.bandwidth = 79; % Fractional bandwidth [%]
        param.radius = 49.57e-3; % Array radius [m]
        param.height = 13.5e-3; % Elevation height [m]
        param.focus = 60e-3; % Elevation focus [m]
    case 'P4-2V'
        % --- P4-2v (Verasonics) ---
        param.fc = 2.72e6; % Transducer center frequency [Hz]
        param.kerf = 50e-6; % Kerf [m]
        param.width = 250e-6; % Width [m]
        param.pitch = 300e-6; % Pitch [m]
        param.Nelements = 64;
        param.bandwidth = 74; % Fractional bandwidth [%]
        param.radius = Inf;
        param.height = 14e-3; % Elevation height [m]
        param.focus = 60e-3; % Elevation focus [m]

% Vermon's matrix arrays (by François Varray)
% note: they contain three void lines in the y-direction
    case 'VERMON-3'
        % --- Vermon 2D probe, 3-MHz ---
        param.fc = 3e6;
        % param.kerf = 50e-6;
        param.width = 250e-6;
        param.height = 250e-6;
        param.bandwidth = 70;
        param.Nelements = 1024;
        %-- Position of the elements (pitch = 300 microns)
        pitch = 300e-6;
        % param.pitch = pitch;
        xe = ((1:32)-16.5)*pitch;
        ye = ([1:8 10:17 19:26 28:35]-18)*pitch;
        [xe,ye] = meshgrid(xe, ye);
        xe = xe'; ye = ye';
        param.elements = [xe(:).'; ye(:).'];
    case 'VERMON-8'
        % --- Vermon 2D probe, 8-MHz ---
        param.fc = 8e6;
        % param.kerf = 50e-6;
        param.width = 250e-6;
        param.height = 250e-6;
        param.bandwidth = 70;
        param.Nelements = 1024;
        %-- Position of the elements (pitch = 300 microns)
        pitch = 300e-6;
        % param.pitch = pitch;
        xe = ((1:32)-16.5)*pitch;
        ye = ([1:8 10:17 19:26 28:35]-18)*pitch;
        [xe,ye] = meshgrid(xe, ye);
        xe = xe'; ye = ye';
        param.elements = [xe(:).'; ye(:).'];

        
       
    %--- From the OLD version of GETPARAM: ---%
    
    case 'PA4-2/20'
        % --- PA4-2/20 ---
        param.fc = 2.5e6; % Transducer center frequency [Hz]
        param.kerf = 50e-6; % Kerf [m]
        param.pitch = 300e-6; % Pitch [m]
        param.height = 14e-3; % Height of element [m]
        param.Nelements = 64;
        param.bandwidth = 60; % Fractional bandwidth [%]
        param.radius = Inf;
    case 'L9-4/38'
        % --- L9-4/38 ---
        param.fc = 5e6; % Transducer center frequency [Hz]
        param.kerf = 35e-6; % Kerf [m]
        param.pitch = 304.8e-6;  % Pitch [m]
        param.height = 6e-3; % Height of element [m]
        param.Nelements = 128;
        param.bandwidth = 65; % Fractional bandwidth [%]
        param.radius = Inf;
    case 'LA530'
        % --- LA530 ---
        param.fc = 3e6; % Transducer center frequency [Hz]
        width = 0.215/1000; %0.2698/1000; % Width of element [m]
        param.kerf = 0.030/1000; %0.035/1000; % Kerf [m]
        param.pitch = width + param.kerf;
        param.height = 6/1000; % Height of element [m]
        param.Nelements = 192;
        param.radius = Inf;
    case 'L14-5/38'
        % --- L14-5/38 ---
        param.fc = 7.2e6; % Transducer center frequency [Hz]
        param.kerf = 25e-6; % Kerf [m]
        param.pitch = 304.8e-6;   % Pitch [m]
        param.height = 4e-3; % Height of element [m]
        param.Nelements = 128;
        param.bandwidth = 70; % Fractional bandwidth [%]
        param.radius = Inf;
    case 'L14-5W/60'
        % --- L14-5W/60 ---
        param.fc = 7.5e6; % Transducer center frequency [Hz]
        param.kerf = 25e-6; % Kerf [m]
        param.pitch = 472e-6;   % Pitch [m]
        param.height = 4e-3; % Height of element [m]
        param.Nelements = 128;
        param.bandwidth = 65; % Fractional bandwidth [%]
        param.radius = Inf;
    case 'P6-3'
        % --- P6-3 ---
        param.fc = 4.5e6; % Transducer center frequency [Hz]
        param.kerf = 25e-6; % Kerf [m]
        param.pitch = 218e-6;   % Pitch [m]
        param.Nelements = 64;
        param.bandwidth = 2/3*100; % Fractional bandwidth [%]
        param.radius = Inf;
        

        
    otherwise
        error(['The probe ' probe ' is unknown.',...
            ' You may complete this function for your own use'])
end


