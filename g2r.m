% This function performs a conversion from ECEF geodetic coordinates to ECEF
% rectangular coordinates.

function v=g2r(c1,lat,c2,lon,h)

v = zeros(3,length(lat));

lat=pi*lat/180;
lon=pi*lon/180;

if c1 == 'S'
    lat=-lat;
elseif c1 ~= 'N'
    error('error, must be N or S');
end

if c2 == 'W'
    lon=-lon;
elseif c2 ~= 'E'
    error('error, must be E or W');
end

% semimajor axis lenght, [m]
a=6378137.0; 

% semiminor axis lenght, [m]
b=6356752.3142;

% flatness of the ellipsoid
f=(a-b)/a;

% eccentricity of the ellipsoid
e=sqrt(f*(2-f));

% length of normal to the ellipsoid, from surface to the 
% ellipsoid to its intersection with the z-axis (see eq. (2.66) in Groves (2008))
N=a./sqrt(1-(e*sin(lat)).^2);

% Calculate the rectangular positions. See eq. (2.70) in Groves (2008).
v(1,:)=(N+h).*cos(lat).*cos(lon); 
v(2,:)=(N+h).*cos(lat).*sin(lon);
v(3,:)=(N.*(1-e^2)+h).*sin(lat);

return