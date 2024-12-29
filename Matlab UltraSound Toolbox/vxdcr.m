function [xV,zV,xe,ze] = vxdcr(param,TXdelay)

%VXDCR   Virtual transducer (XDCR)
%   [xV,zV] = VXDCR(PARAM,DELAYS) returns the positions of the elements
%   of the virtual transducer.
%
%   [xV,zV,xE,zE] = VXDCR(...) also returns the positions of the elements
%   of the actual transducer.
%
%   The virtual transducer is an "equivalent" transducer assuming null
%   delays.
%
%
%   Examples:
%   --------
%   % Example #1: Phased-array transducer with focused transmission
%   figure
%   param = getparam('P4-2v');
%   dels = txdelay(0,3e-2,param);
%   [xv,zv,xe,ze] = vxdcr(param,dels);
%   plot(xv*1e3,zv*1e3,'ro',xe*1e3,ze*1e3,'bs')
%   axis equal ij tight
%   xlabel('[mm]')
%   legend({'virtual transducer','actual transducer'},'Location','southoutside')
%   title('64-element phased array w/ focused transmission')
%
%   % Example #2: Convex transducer with plane-wave transmission
%   figure
%   param = getparam('C5-2v');
%   tilt = 10/180*pi; % tilt of 10 degrees
%   dels = txdelay(param,tilt);
%   [xv,zv,xe,ze] = vxdcr(param,dels);
%   plot(xv*1e3,zv*1e3,'ro',xe*1e3,ze*1e3,'bs')
%   axis equal ij tight
%   xlabel('[mm]')
%   legend({'virtual transducer','actual transducer'},'Location','southoutside')
%   title('128-element convex array w/ plane-wave transmission')
%
%   % Example #3: Convex transducer with a sub-aperture & focused wave
%   figure
%   param = getparam('C5-2v');
%   dels = NaN(1,param.Nelements);
%   param_sub = param; % we'll use a sub-aperture of 24 elements
%   param_sub.Nelements = 21;
%   dels(11:31) = txdelay(0,1.5e-2,param_sub);
%   [xv,zv,xe,ze] = vxdcr(param,dels);
%   plot(xv*1e3,zv*1e3,'ro',xe*1e3,ze*1e3,'bs')
%   axis equal ij tight
%   xlabel('[mm]')
%   legend({'virtual transducer','actual transducer'},'Location','southoutside')
%   title('128-element convex array w/ sub-aperture focused transmission')
%
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also GETPARAM, TXDELAY, VIEWXDCR, DASMTX.
%
%   -- Damien Garcia -- 2023/04, last update: 2024/05/29
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>

narginchk(1,2)
assert(isstruct(param),'PARAM must be a structure.')
% param = IgnoreCaseInFieldNames(param);
if nargin==1
    assert(isfield(param,'TXdelay'),...
        'A TX delay vector (PARAM.TXdelay or DELAYS) is required.')
    TXdelay = param.TXdelay;
elseif isfield(param,'TXdelay')
    assert(isequal(TXdelay,param.TXdelay),...
        'If both specified, PARAM.TXdelay and DELAYS must be equal.')
end
assert(all(TXdelay(~isnan(TXdelay))>=0),'The delays must be non-negative.')

%-- Number of elements
nc = size(TXdelay,2);
% Note: param.Nelements can be required in other functions of the
%       Matlab Ultrasound Toolbox
if isfield(param,'Nelements')
    assert(param.Nelements==nc,...
        'PARAM.TXdelay or DELAYS must be of length PARAM.Nelements.')
end

%-- Longitudinal velocity (in m/s)
if ~isfield(param,'c')
    param.c = 1540;
end
c = param.c;

%-- Centers of the tranducer elements (x- and z-coordinates)
if isinf(param.radius)
    isLINEAR = true;
    % Linear array
    xe = ((0:nc-1)-(nc-1)/2)*param.pitch;
    ze = zeros(size(xe));
else
    isLINEAR = false;
    % Convex array
    R = param.radius;
    % https://en.wikipedia.org/wiki/Circular_segment
    chord = 2*R*sin(asin(param.pitch/2/R)*(nc-1));
    h = sqrt(R^2-chord^2/4);
    THe = linspace(atan2(-chord/2,h),atan2(chord/2,h),nc);
    % THe = angle of the normal to element #e about the y-axis
    ze = R*cos(THe);
    xe = R*sin(THe);
    ze = ze-h;
    % Note: the center of the circular segment is then (0,-h)
end

%-- First check which elements were transmitting
WasTransmitting = ~isnan(TXdelay);
assert(sum(abs(diff(WasTransmitting)))<3,...
    'Multiple transmitting sub-apertures are not allowed during beamforming; separate them accordingly.')
nTX = nnz(WasTransmitting); % number of transmitting elements
assert(nTX>=3,...
    'The number of neighboring transmitting elements must be at least 3.')

xeT = xe(WasTransmitting);
zeT = ze(WasTransmitting);
TXdelay = TXdelay(WasTransmitting);

T2 = TXdelay.^2;
dT2 = diff2(xeT,T2);
err = min(diff(xeT))*1e-3;

%-- Virtual transducer
if isLINEAR
    xV = xeT - 0.5*c^2*dT2;
    zV = -sqrt(abs(c^2*T2-(xV-xeT).^2));
else
    % for a convex array
    xV = xeT;
    dzedxe = diff2(xeT,zeT);
    tmp = sqrt(abs(c^2*T2-(xV-xeT).^2));
    opt = optimset('TolFun',err,'TolX',err);
    for k = 1:length(xV)
        myfun = @(xv) (xeT(k)-xv) +...
            tmp(k).*dzedxe(k) +...
            -0.5*c^2*dT2(k);
        xV(k) = fzero(@(xv) myfun(xv),xV(k),opt);
    end
    zV = zeT-tmp;
end

%-- Inappropriate delays can generate incorrect virtual transducers!
test = all(diff(xV)>0) & all((c^2*TXdelay.^2-(xV-xeT).^2)>-err);
if ~test
    warning('off','backtrace')
    warning('The delays do not allow a virtual transducer to be calculated.')
    xV = NaN(size(xV)); zV = NaN(size(zV));
end

sgn_d2zV = sign(diff(diff(zV,2)));
test = all(ismember(sgn_d2zV,[0 1])) | all(ismember(sgn_d2zV,[0 -1]));
if ~test
    warning('off','backtrace')
    warning('The virtual transducer is not convex or concave.')
end



end



function dy = diff2(x,y)
% 1st derivative of Y using a 2nd order difference method
% [simplified version: for 1-D arrays]

m = length(y);
dx = diff(x);
dy = zeros(size(y));

% y'(x1)
dy(1) = (1/dx(1)+1/dx(2))*(y(2)-y(1))+...
    dx(1)/(dx(1)*dx(2)+dx(2)^2)*(y(1)-y(3));
% y'(xm)
dy(m) = (1/dx(m-2)+1/dx(m-1))*(y(m)-y(m-1))+...
    dx(m-1)/(dx(m-1)*dx(m-2)+dx(m-2)^2)*(y(m-2)-y(m));
% y'(xi) (i>1 & i<m)
dx1 = dx(1:m-2);
dx2 = dx(2:m-1);
y1 = y(1:m-2); y2 = y(2:m-1); y3 = y(3:m);
dy(2:m-1) = 1./(dx1.*dx2.*(dx1+dx2)).*...
    (-dx2.^2.*y1+(dx2.^2-dx1.^2).*y2+dx1.^2.*y3);

end


