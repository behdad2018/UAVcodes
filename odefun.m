
function dydt = odefun(t,y)

dydt = zeros(3,1);

M

L

dydt = -inv(M)*L + f
dydt(1) = y(1)+2*y(2);
dydt(2) = 3*y(1)+2*y(2);

end