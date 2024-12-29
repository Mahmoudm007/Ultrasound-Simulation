function p = viewxdcr(param,unit)

%VIEWXDCR   View your transducer (XDCR)
%   VIEWXDCR(PARAM) plots the transducer (xdcr) defined by the parameters
%   included in the structure PARAM.
%
%   The transducer can be a rectilinear or convex ULA (uniform linear
%   array), for a use with TXDELAY, PFIELD and SIMUS, or a planar array
%   such as a matrix array, for a use with TXDELAY3, PFIELD3 and SIMUS3.
%   See below the fields of the PARAM structure.
%
%   VIEWXDCR(PARAM,UNIT) uses UNIT as a unit distance. UNIT can be 'm',
%   'cm', or 'mm' (it is the default). UNIT can also be scalar S, in which
%   case the unit is S*m.
%
%   P = VIEWXDCR(PARAM) also returns the patch object that contains the
%   data for all the elements of the transducer.
%
%   Note: If a transmit apodization is given (PARAM.TXapodization), it
%   determines the element colors.
%
%   PARAM is a structure that contains the following fields:
%   -------------------------------------------------------
%       *** UNIFORM LINEAR ARRAY (rectilinear or convex) ***
%   1) PARAM.Nelements: number of elements in the transducer array (REQUIRED)
%   2) PARAM.pitch: pitch of the array (in m, REQUIRED)
%   3) PARAM.width: element width (in m, REQUIRED)
%       or PARAM.kerf: kerf width (in m, REQUIRED)
%   4) PARAM.height: element height (in m, default = Inf: 2-D acoustics) 
%   5) PARAM.radius: radius of curvature (in m, default = Inf)
%
%       *** 2-D PLANAR ARRAY (e.g. matrix array) ***
%   6) PARAM.elements: x- and y-coordinates of the element centers
%       (in m, REQUIRED). It MUST be a two-row matrix, with the 1st and 2nd
%       rows containing the x and y coordinates, respectively.
%   7) PARAM.width: element width, in the x-direction (in m, REQUIRED)
%   8) PARAM.height: element height, in the y-direction (in m, REQUIRED)
%
%       *** TRANSMIT PARAMETERS ***
%   9)  PARAM.TXapodization: transmit apodization (default: no apodization)
%
%
%   Examples:
%   --------
%   % Phased-array transducer with apodization
%   figure
%   param = getparam('P4-2v');
%   param.TXapodization = sin(pi*linspace(0,63,64)/64).^2;
%   viewxdcr(param)
%   colormap parula
%   colorbar
%   title('64-element phased array w/ apodization')
%
%   % Convex transducer
%   figure
%   param = getparam('C5-2v');
%   viewxdcr(param)
%   title('128-element convex array')
%
%   % Matrix array with 32x32 elements with apodization
%   figure
%   param = [];
%   param.width = 250e-6;
%   param.height = 250e-6;
%   % position of the elements (pitch = 300 microns)
%   pitch = 300e-6;
%   [xe,ye] = meshgrid(((1:32)-16.5)*pitch);
%   param.elements = [xe(:).'; ye(:).'];
%   % transmit apodization
%   apod = sin(pi*linspace(0,31,32)/32).^2;
%   apod = apod'*apod;
%   param.TXapodization = apod(:);
%   % view the transducer
%   viewxdcr(param)
%   title('32\times32-element matrix array w/ apodization')
%
%
%   This function is part of <a
%   href="matlab:web('https://www.biomecardio.com/MUST')">MUST</a> (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   See also GETPARAM, PFIELD, PFIELD3, SIMUS, SIMUS3.
%
%   -- Damien Garcia -- 2022/11, last update: 2023/01/20
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>


narginchk(1,2)
assert(isstruct(param),'PARAM must be a structure.')
if nargin==1, unit = 'mm'; end
assert(isscalar({unit}),'Wrong argument for UNIT.')
if ~isnumeric(unit)
    assert(ismember(unit,{'m','cm','mm'}),...
        'Wrong argument for UNIT.')
    zlabelname = unit;
    switch unit
        case 'm'
            unit = 1;
        case 'cm'
            unit = 100;
        case 'mm'
            unit = 1000;
    end
elseif isscalar(unit) && isnumeric(unit)
    assert(unit>0,'If a scalar, UNIT must be >0.')
    zlabelname = [num2str(unit) '\timesm'];
    unit = 1/unit;
end

isULA = ~isfield(param,'elements');

%-- Element height (in m)
if ~isfield(param,'height')
    param.height = Inf; % default = line array
end
ElementHeight = param.height;
assert(isnumeric(ElementHeight) && isscalar(ElementHeight) &&...
    ElementHeight>0,'The element height must be positive.')

if isULA % we have a Uniform Linear Array (rectilinear or convex)

    %-- Number of elements
    assert(isfield(param,'Nelements'),'The number of elements (PARAM.Nelements) is required.')
    NumberOfElements = param.Nelements;
    assert(isnumeric(NumberOfElements) && isscalar(NumberOfElements) &&...
        NumberOfElements==abs(round(NumberOfElements)),...
        'The number of elements is not conventional.')

    %-- Pitch (in m)
    assert(isfield(param,'pitch'),'A pitch value (PARAM.pitch) is required.')
    pitch = param.pitch;
    assert(isnumeric(pitch) && isscalar(pitch) && pitch>0,...
        'The pitch must be positive.')

    %-- Element width and/or Kerf width (in m)
    if isfield(param,'width') && isfield(param,'kerf')
        assert(abs(pitch-param.width-param.kerf)<eps,...
            'The pitch must be equal to (kerf width + element width).')
    elseif isfield(param,'kerf')
        assert(isnumeric(param.kerf) && isscalar(param.kerf) &&...
            param.kerf>0,'The kerf width must be positive.')
        param.width = pitch-param.kerf;
    elseif isfield(param,'width')
        assert(isnumeric(param.width) && isscalar(param.width) &&...
            param.width>0,'The element width must be positive.')
        param.kerf = pitch-param.width;
    else
        error(['An element width (PARAM.width) ',...
            'or kerf width (PARAM.kerf) is required.'])
    end
    ElementWidth = param.width;

    %-- Radius of curvature (in m) - convex array
    if ~isfield(param,'radius')
        param.radius = Inf; % default = linear array
    end
    RadiusOfCurvature = param.radius;
    assert(isnumeric(RadiusOfCurvature) && isscalar(RadiusOfCurvature) &&...
        RadiusOfCurvature>0,'The radius of curvature must be positive.')

    %-- Concavity
    if isfield(param,'concave')
        assert(islogical(param.concave),'PARAM.concave must be logical.')
    else
        param.concave = false;
    end

else % we have a 2-D planar array

    %-- Transducer elements
    assert(isfield(param,'elements'),...
        ['PARAM.elements must contain the x- and y-locations ',...
        'of the transducer elements.'])
    assert(size(param.elements,1)==2,...
        ['PARAM.elements must have two rows that contain the x (1st row) ',...
        'and y (2nd row) coordinates of the transducer elements.'])
    
    %-- Element width
    if isfield(param,'width')
        assert(isnumeric(param.width) && isscalar(param.width) &&...
            param.width>0,'The element width must be positive.')
    else
        error('An element width (PARAM.width) is required.')
    end
    ElementWidth = param.width;

    NumberOfElements = size(param.elements,2);
    param.concave = false;

end

%-- Coordinates of the element centers
if isULA % this is a uniform linear array (rectilinear or convex)

    %-- Centers of the tranducer elements (x- and z-coordinates)
    if isinf(RadiusOfCurvature)
        % we have a LINEAR ARRAY
        xe = ((0:NumberOfElements-1)-(NumberOfElements-1)/2)*param.pitch;
        ze = zeros(1,NumberOfElements);
        THe = zeros(1,NumberOfElements);
    else
        % we have a CONVEX ARRAY
        chord = 2*RadiusOfCurvature*...
            sin(asin(param.pitch/2/RadiusOfCurvature)*(NumberOfElements-1));
        h = sqrt(RadiusOfCurvature^2-chord^2/4); % apothem
        % https://en.wikipedia.org/wiki/Circular_segment
        % THe = angle of the normal to element #e about the y-axis
        THe = linspace(atan2(-chord/2,h),atan2(chord/2,h),NumberOfElements);
        ze = RadiusOfCurvature*cos(THe);
        xe = RadiusOfCurvature*sin(THe);
        ze = (ze-h) + param.concave*(RadiusOfCurvature+h-2*ze);
        % Note: the center of the circular segment is then (0,-h)
    end
    ye = zeros(1,numel(xe));

else % this is a 2-D planar array
    xe = param.elements(1,:);
    ye = param.elements(2,:);
    ze = zeros(size(xe));
    THe = zeros(size(xe));
end

%-- Coordinates of the rectangles for the function PATCH
if isinf(ElementHeight), ElementHeight=0; end
xp = xe + cos(THe).*[-1 1 1 -1]'*ElementWidth/2;
yp = ye + [-1 -1 1 1]'*ElementHeight/2;
zp = ze + (1-2*param.concave)*sin(THe).*[1 -1 -1 1]'*ElementWidth/2;

%-- Colors of the polygons (PARAM.TXapodization)
if ~isfield(param,'TXapodization')
    C = ones(1,NumberOfElements);
else
    assert(isvector(param.TXapodization) && isnumeric(param.TXapodization),...
        'PARAM.TXapodization must be a vector')
    assert(numel(param.TXapodization)==NumberOfElements,...
        'PARAM.TXapodization must be of length = (number of elements)')
    C = param.TXapodization;
end

%-- Patches
if nargout>0
    p = patch(xp*unit,yp*unit,zp*unit,C,'EdgeColor','k');
else
    patch(xp*unit,yp*unit,zp*unit,C,'EdgeColor','k')
end
view(3)
set(gca,'DataAspectRatio',[1 1 1],'zdir','reverse',...
    'box','on','boxstyle','full')
axis equal
zlim([-1e-3 1e-3]*unit+[min(ze) max(ze)]*unit)
zlabel(['[ ' zlabelname ' ]'])


