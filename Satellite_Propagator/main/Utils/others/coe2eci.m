function eci_el = coe2eci(coe_el,mu)
% Algorithm 10 from Vallado
a = coe_el(1);
e = coe_el(2);
i = coe_el(3);
w = coe_el(4);
ohm = coe_el(5);
nu = coe_el(6);

% mu = 398600.4418;
% TODO: The agorithm is not completely implemented
% Step 1:
% Calculate the semi perimeter p
if e ==1
    error('e is 1');
end
p = a*(1-e^2);
r_PQW = [p*cos(nu)/(1+e*cos(nu)); p*sin(nu)/(1+e*cos(nu)); 0];
v_PQW = [-sqrt(mu/p)*sin(nu); sqrt(mu/p)*(e+cos(nu)); 0];
temp_mat = [cos(ohm)*cos(w)-sin(ohm)*sin(w)*cos(i), -cos(ohm)*sin(w)-sin(ohm)*cos(w)*cos(i), sin(ohm)*sin(i);...
            sin(ohm)*cos(w)+cos(ohm)*sin(w)*cos(i), -sin(ohm)*sin(w)+cos(ohm)*cos(w)*cos(i), -cos(ohm)*sin(i);...
            sin(w)*sin(i), cos(w)*sin(i), cos(i)];
r_IJK = temp_mat*r_PQW;
v_IJK = temp_mat*v_PQW;
        
eci_el = [r_IJK',v_IJK'];

end