% This function calculates the gravity force.

function g=gravity(lat)

lat=pi/180*lat;

e = 0.0818191908425;                                                % See page 38 in Groves (2008).
g=9.7803253359*(1+0.001931853*sin(lat)^2)/sqrt(1-e^2*sin(2*lat)^2); % See eq. (2.85) in Groves (2008).

return;




