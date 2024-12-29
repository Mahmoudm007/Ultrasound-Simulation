function varargout = rmscat(varargin)

%RMSCAT   Remove insignificant scatterers for faster simulations
%   [X2,Y2,Z2,RC2] = RMSCAT(X,Y,Z,RC,DELAYS,PARAM) removes the scatterers
%   that would have little effect on an ultrasound image simulated by SIMUS
%   or SIMUS3 and of dynamic range 40 dB. The retained scatterers are
%   located at {X2,Y2,Z2} and their reflection coefficients are RC2.
%
%   RMSCAT can be used before SIMUS or SIMUS3. After using RMSCAT, you can
%   thus consider running RF = SIMUS(X2,Y2,Z2,RC2,DELAYS,PARAM) instead of
%   RF = SIMUS(X,Y,Z,RC,DELAYS,PARAM). Same for SIMUS3.
%   
%   [X2,Z2,RC2] = RMSCAT(X,Z,RC,DELAYS,PARAM) assumes that you will use
%   2-D acoustics in the X-Z plane: RF = SIMUS(X2,Z2,RC2,DELAYS,PARAM) or
%   RF = SIMUS(X2,[],Z2,RC2,DELAYS,PARAM).
%
%   I = RMSCAT(...) returns the linear indices corresponding to the
%   retained scatterers, i.e. X2 = X(I), Y2 = Y(I),...
%
%   [...] = RMSCAT(...,DR) removes the scatterers that would have little
%   effect on a simulated ultrasound of dynamic range DR dB.
%   By default, DR = 40.
%
%
%   Example: cardiac scatterers:
%   ---------------------------
%   %-- Read the heart image and make it gray
%   I = imread('heart.jpg');
%   I = rgb2gray(I);
%   %-- Parameters of the cardiac phased array
%   param = getparam('P4-2v');
%   %-- Pseudorandom distribution of scatterers (depth is 15 cm)
%   [xs,~,zs,RC] = genscat([NaN 15e-2],param,I);
%   %-- Display the scatterers in a dB scale
%   subplot(131)
%   scatter(xs*1e2,zs*1e2,2,20*log10(RC/max(RC(:))),'filled')
%   caxis([-40 0])
%   colormap(dopplermap)
%   axis equal ij tight
%   V = axis;
%   set(gca,'XColor','none','box','off')
%   title([int2str(numel(RC)) ' cardiac scatterers'])
%   ylabel('[cm]')
%   %-- RMS pressure field for a focused wave with a cardiac phased array
%   subplot(132)
%   param = getparam('P4-2v');
%   dels = txdelay(-1e-2,8e-2,param); % focus at (-1,8) cm
%   [x,z] = meshgrid(linspace(-7.5,7.5,256)*1e-2,linspace(.1,15,256)*1e-2);
%   P = pfield(x,0*x,z,dels,param);
%   imagesc(x(1,:)*1e2,z(:,1)*1e2,20*log10(P/max(P(:))))
%   axis equal ij tight
%   caxis([-40 0])
%   title('acoustic field ([-40,0] dB)')
%   %-- Retain the significant scatterers this sequence
%   subplot(133)
%   [xs2,~,zs2,RC2] = rmscat(xs,0*xs,zs,RC,dels,param); 
%   scatter(xs2*1e2,zs2*1e2,2,20*log10(RC2/max(RC2(:))),'filled')
%   caxis([-40 0])
%   colormap(dopplermap)
%   axis equal ij, axis(V)
%   set(gca,'XColor','none','box','off')
%   title([int2str(numel(RC2)) ' cardiac scatterers'])
%   
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also GENSCAT, SIMUS, SIMUS3, CITE.
%
%   -- Damien Garcia -- 2024/02, last update 2024/02/22
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>

narginchk(5,7)
nargoutchk(0,5)
if nargout==2
    error('The number of output arguments must be 1, 3, or 4.')
end

if isstruct(varargin{nargin})
    % rmscat(...,param)
    NArg = nargin-1;
    param = varargin{end};
    DR = 40;
elseif isstruct(varargin{nargin-1})
    % rmscat(...,param,DR)
    NArg = nargin-2;
    param = varargin{end-1};
    DR = varargin{end};
else
    error('The PARAM structure is missing.')
end
assert(DR>0,'The dynamic range DR (in dB) must be >0')

xs = varargin{1};
switch NArg
    case 4
        % rmscat(xs,zs,RC,TXdelays,...)
        ys = [];
        zs = varargin{2};
        RC = varargin{3};
        TXdelays = varargin{4};
    case 5
        % rmscat(xs,ys,zs,RC,TXdelays,...)
        ys = varargin{2};
        zs = varargin{3};
        RC = varargin{4};
        TXdelays = varargin{5};
end

opt.WaitBar = false;

% Adjust options and parameters to make PFIELD or PFIELD3 faster
% (no need for high precision)
opt.dBThresh = -20;
param.bandwidth = 20;
param.TXnow = 4;
opt.ElementSplitting = 1;

if isfield(param,'elements')
    P = pfield3(xs,ys,zs,TXdelays,param,opt);
else
    P = pfield(xs,ys,zs,TXdelays,param,opt);
end

PRC = P.*RC;
clear P
idx = PRC>max(PRC,[],'all')*10^(-DR/20);

if nargout==1
    varargout{1} = idx;
elseif nargout==3
    varargout{1} = xs(idx);
    varargout{2} = zs(idx);
    varargout{3} = RC(idx);
else
    varargout{1} = xs(idx);
    if isempty(ys)
        varargout{2} = [];
    else
        varargout{2} = ys(idx);
    end
    varargout{3} = zs(idx);
    varargout{4} = RC(idx);
end