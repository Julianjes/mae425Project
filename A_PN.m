function Apn = A_PN(COE)
% Input: Structure with classical orbital elements

omega = COE.RAAN;
w = COE.w;
i = COE.i;

Apn= [cosd(omega)*cosd(w)-sind(omega)*sind(w)*cosd(i)  -cosd(omega)*sind(w)-sind(omega)*cosd(w)*cosd(i)  sind(omega)*sind(i);
       sind(omega)*cosd(w)+cosd(omega)*sind(w)*cosd(i)  -sind(omega)*sind(w)+cosd(omega)*cosd(w)*cosd(i)  -cosd(omega)*sind(i);
                       sind(w)*sind(i)                           cosd(w)*sind(i)                                cosd(i)];


end