function r_ECEF = ECI_ECEF(JD,r_ECI)
% JD is Julian Date
% r_ECI is position in ECI Frame

theta_ERA = mod(286.46061837504 + 360.985612288808*(JD - 2451545.0),360);

rot3 = [cosd(theta_ERA) sind(theta_ERA) 0
        -sind(theta_ERA) cosd(theta_ERA) 0
            0               0              1];
r_ECEF = rot3*r_ECI;

end