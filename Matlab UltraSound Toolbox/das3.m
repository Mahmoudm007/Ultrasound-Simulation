function [bfSIG,param] = das3(SIG,x,y,z,varargin)

%DAS3   Delay-And-Sum beamforming of RF or I/Q volume signals
%   BFSIG = DAS3(SIG,X,Y,Z,DELAYS,PARAM) beamforms the RF or I/Q volume
%   signals stored in the array SIG, and returns the beamformed signals
%   BFSIG. The signals are beamformed at the points specified by X, Y, and
%   Z.
%
%   >--- Try it: enter "das3" in the command window for an example ---<
%
%     1) SIG must be a 2-D array. The first dimension (i.e. each
%        column) corresponds to a single RF or I/Q signal over (fast-)
%        time, with the COLUMN #n corresponding to the ELEMENT #n of
%        PARAM.elements.
%     2) DELAYS are the transmit time delays (in s). One must have
%        numel(DELAYS) = size(SIG,2). If a sub-aperture was used during
%        transmission, use DELAYS(i) = NaN if element #i of the array was
%        off.
%     3) PARAM is a structure that contains the parameter values required
%        for the delay-and-sum (see below for details).
%
%   Note: SIG must be complex when DAS-beamforming I/Q data
%         (i.e. SIG = complex(I,Q) = I + 1i*Q).
%
%   DAS3(SIG,X,Y,Z,PARAM) uses DELAYS = PARAM.TXdelay.
%
%   DAS(...,METHOD) specifies the interpolation method. The available
%   methods are decribed in NOTE #3 below.
%
%   [BFSIG,PARAM] = DAS3(...) also returns the structure PARAM.
%
%   ---
%   NOTE #1: DAS3 uses a standard diffraction summation (delay-and-sum). It
%   calls the function DASMTX3.
%   ---
%   NOTE #2: Interpolation method
%   By default DAS3 uses a linear interpolation to generate the DAS matrix.
%   To specify the interpolation method, use DAS3(...,METHOD), with METHOD
%   being:
%      'nearest'   - nearest neighbor interpolation
%      'linear'    - (default) linear interpolation
%      'quadratic' - quadratic interpolation
%      'lanczos3'  - 3-lobe Lanczos interpolation
%      '5points'   - 5-point least-squares parabolic interpolation
%      'lanczos5'  - 5-lobe Lanczos interpolation
%
%   The linear interpolation (it is a 2-point method) returns a matrix
%   twice denser than the nearest-neighbor interpolation. It is 3, 4, 5, 6
%   times denser for 'quadratic', 'lanczos3', '5points', 'lanczos5',
%   respectively (they are 3-to-6-point methods).
%   ---
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%   1)  PARAM.fs: sampling frequency (in Hz, REQUIRED)
%   2)  PARAM.elements: coordinates of the transducer elements (in m, REQUIRED)
%       PARAM.elements must contain the x- and y-coordinates of the
%       transducer elements (the z-coordinates are zero if not given). It
%       must be a matrix with 2 (or 3) rows corresponding to the
%       x- and y-coordinates (and optionally z), respectively.
%   3)  PARAM.fc: center frequency (in Hz, REQUIRED for I/Q signals)
%   4)  PARAM.TXdelay: transmission delays (in s, required if DELAYS is not given)
%   5)  PARAM.c: longitudinal velocity (in m/s, default = 1540 m/s)
%   6)  PARAM.t0: start time for reception (in s, default = 0 s)
%   7)  PARAM.fnumber: reception f-number (default = [0 0], i.e. full aperture)
%       PARAM.fnumber(1) = reception f-number in the azimuthal x-direction.
%       PARAM.fnumber(2) = reception f-number in the elevation y-direction.
%
%   Passive imaging
%   ---------------
%   8)  PARAM.passive: must be true for passive imaging (i.e. no transmit).
%       The default is false.
%
%   If you need to beamform a large series of ultrasound signals acquired
%   with a same probe and a same transmit sequence, DASMTX3 is recommended.
%   BEWARE, however. DASMTX3 can generate tall sparse DAS matrices when
%   beamforming large volume data! Consider chunking your datasets.
%
%   DAS3 chunks the input data (X,Y,Z) to avoid using too large matrices.
%
%
%   REFERENCES:
%   ----------
%   V Perrot, M Polichetti, F Varray, D Garcia. So you think you can DAS? A
%   viewpoint on delay-and-sum beamforming. Ultrasonics 111, 106309. <a
%   href="matlab:web('https://www.biomecardio.com/publis/ultrasonics21.pdf')">PDF here</a>
%
%
%   Example:
%   -------
%   %-- Plane wave w/ a matrix-array and PSF at 3 cm depth 
%   %- 3-MHz matrix array with 32x32 elements
%   param = [];
%   param.fc = 3e6;
%   param.bandwidth = 70;
%   param.width = 250e-6;
%   param.height = 250e-6;
%   %- Position of the elements (pitch = 300 microns)
%   pitch = 300e-6;
%   [xe,ye] = meshgrid(((1:32)-16.5)*pitch);
%   param.elements = [xe(:).'; ye(:).'];
%   %- PSF position
%   x0 = 0; y0 = 0; z0 = 3e-2;
%   %- Transmit time delays (plane wave)
%   dels = zeros(1,1024);
%   %- Simulate RF signals with SIMUS3
%   [RF,param] = simus3(x0,y0,z0,1,dels,param);
%   %- Demodulate
%   IQ = rf2iq(RF,param);
%   %- 3-D grid for beamforming
%   lambda = 1540/param.fc;
%   [xi,yi,zi] = meshgrid(-2e-2:lambda:2e-2,-2e-2:lambda:2e-2,...
%       2.5e-2:lambda/2:3.1e-2);
%   %-- Beamform using DAS3
%   IQb = das3(IQ,xi,yi,zi,dels,param);
%   %-- Log-envelope
%   I = 20*log10(abs(IQb)/max(abs(IQb(:))));
%   %- Figure of the PSF
%   figure
%   I(1:round(size(I,1)/2),1:round(size(I,2)/2),:) = NaN;
%   for k = [-40:10:-10 -5 -1]
%       isosurface(xi*1e2,yi*1e2,zi*1e2,I,k)
%   end
%   colormap([1-hot;hot])
%   colorbar
%   box on, grid on
%   zlabel('[cm]')
%   title('PSF at (0,0,3) cm [dB]')
%
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also DASMTX3, DAS, TXDELAY3, RF2IQ, SIMUS3.
%
%   -- Damien Garcia -- 2022/12, last update 2023/03/01
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>


if nargin==0
    if nargout>0
        [bfSIG,param] = RunTheExample;
    else
        RunTheExample;
    end
    return
end

assert(nargin>4,'Not enough input arguments.')
assert(nargin<8,'Too many input arguments.')

siz0 = size(x);
assert(isequal(siz0,size(y),size(z)),'X, Y, and Z must be of same size.')
assert(ndims(SIG)<=2,['SIG must be a 2D array whose each column ',...
    'corresponds to an RF or IQ signal acquired by a single element']);

[nl,nc] = size(SIG);

%-- check if we have I/Q signals
isIQ = ~isreal(SIG);

%-- Check METHOD
if ischar(varargin{end})
    method = varargin{end};
else
    method = 'linear';
end
tmp = strcmpi(method,...
    {'nearest','linear','quadratic','lanczos3','5points','lanczos5'});
if ~any(tmp)
    error(['METHOD must be ''nearest'', ''linear'', ''quadratic'',',...
        ' ''Lanczos3'', ''5points'' or ''Lanczos5''.'])
end
Npoints = find(tmp);

%-- Check for potential errors
try
    dasmtx3([nl nc],x(1),y(1),z(1),varargin{:});
catch ME
    throw(ME)
end



%---- Chunking ----
%
% Large volume data can generate tall DAS matrices. The input data (X,Y,Z)
% are chunked to avoid out-of-memory issues.

%-- Maximum possible array bytes
if ispc
    % The memory function is available only on Microsoft Windows platforms
    mem = memory;
    MPAB = mem.MaxPossibleArrayBytes;
end

% The number of bytes required to store a sparse matrix M is roughly:
%     bytes = 16*nnz(M) + 8*(size(M,2)+1)
% (for a 64-bit system)
%
% In our case:
% nnz(M) < (number of transducer elements)*...
%          (number of interpolating points)*...
%          (number of grid points)
% size(M,2) = min(number of RF/IQ samples,number of grid points)
% 
% Roughly:
% bytes < 16*(number of transducer elements)*...
%            (number of interpolating points)*...
%            (number of grid points)
%

%-- Number of chunks
NoE = size(SIG,2); % number of elements
bytes = 16*NoE*Npoints*numel(x);
factor = 20; % other large variables in DASMTX3 + ...
             % compromise for-loops vs. memory
             % (arbitrary value > 4-5)
if ispc
    Nchunks = ceil(factor*bytes/MPAB);
else
    Nchunks = 10; % arbitrary value
end

%---- end of Chunking ----



SIG = SIG(:);
bfSIG = zeros(siz0,'like',SIG);
idx = round(linspace(1,numel(x)+1,Nchunks+1));

for k = 1:Nchunks
    %-- DAS matrices using DASMTX3
    [M,param] = dasmtx3((~isIQ*1+isIQ*1i)*[nl nc],...
        x(idx(k):idx(k+1)-1),...
        y(idx(k):idx(k+1)-1),...
        z(idx(k):idx(k+1)-1),...
        varargin{:});

    %-- Delay-and-Sum
    bfSIG(idx(k):idx(k+1)-1) = M*SIG;
end

bfSIG = reshape(bfSIG,siz0);

end


function [IQb,param] = RunTheExample

%-- 3-MHz matrix array with 32x32 elements
param = [];
param.fc = 3e6;
param.bandwidth = 70;
param.width = 250e-6;
param.height = 250e-6;

%-- Position of the elements (pitch = 300 microns)
pitch = 300e-6;
[xe,ye] = meshgrid(((1:32)-16.5)*pitch);
param.elements = [xe(:).'; ye(:).'];

%-- PSF position
x0 = 0; y0 = 0; z0 = 3e-2;

%-- Transmit time delays using TXDELAY3
% dels = txdelay3(x0,y0,z0,param);
dels = zeros(1,1024);

%-- Pressure field with PFIELD3
n = 24;
[xi,yi,zi] = meshgrid(linspace(-5e-3,5e-3,n),linspace(-5e-3,5e-3,n),...
    linspace(0,6e-2,4*n));
RP = pfield3(xi,yi,zi,dels,param);

% Figure of the pressure field
figure
%--
% y-z slice @ x = 0 cm
% x-z slice @ y = 0 cm
% x-y slice @ z = 3 cm
slice(xi*1e2,yi*1e2,zi*1e2,RP,0,0,3)
%--
set(gca,'zdir','reverse')
zlabel('[mm]')
shading flat
axis equal
colormap([1-hot;hot])
hold on
plot3(xe(:)*1e2,ye(:)*1e2,0*xe(:),'.')
title('A radiation pattern with a 32{\times}32 matrix array')

%-- Simulate RF signals with SIMUS3
[RF,param] = simus3(x0,y0,z0,1,dels,param);

%-- Demodulate
IQ = rf2iq(RF,param);

%-- 3-D grid for beamforming
lambda = 1540/param.fc;
[xi,yi,zi] = meshgrid(-2e-2:lambda:2e-2,-2e-2:lambda:2e-2,...
    2.5e-2:lambda/2:3.2e-2);

%-- Beamform
IQb = das3(IQ,xi,yi,zi,dels,param);

%-- Log-envelope
I = 20*log10(abs(IQb)/max(abs(IQb(:))));

%-- Figure of the PSF
figure
I(1:round(size(I,1)/2),1:round(size(I,2)/2),:) = NaN;
for k = [-40:10:-10 -5 -1]
    isosurface(xi*1e2,yi*1e2,zi*1e2,I,k)
end
colormap([1-hot;hot])
colorbar
box on, grid on
zlabel('[cm]')
title('PSF at (0,0,3) cm [dB]')

end


