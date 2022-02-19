% This function calculates the rotation matrix for rotation from e frame to t frame.

function C=Ce2t(c1,lat,c2,lon)

if(c1=='S')
    lat=-lat;
end

if(c2=='W')
    lon=-lon;
end

lat=lat/180*pi;
lon=lon/180*pi;

clat=cos(lat);
slat=sin(lat);

clon=cos(lon);
slon=sin(lon);

% See eq. (2.99) in Groves (2008).
C=[-slat*clon -slat*slon  clat; 
   -slon          clon      0 ;
   -clat*clon -clat*slon -slat];

end
