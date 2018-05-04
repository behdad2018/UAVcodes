% performace propeller model, Written by Behdad Davoudi, Computaional
% Aeroscience Lab, Aerospace Engineering Department, University of
% Michigan, 2017


% 
% rpm=5000;
% a=2*pi;
% R=0.2;
% c=0.05;
% ro=1;
% Nb=2;
% A=pi*R^2;
% sig=Nb*c/(pi*R);
% om=rpm*(2*pi/60);
% gam=ro*a*c*R^4/1;
% Vt=om*R;
% th0=degtorad(20);
% alpha=degtorad(90); % for a propeller it is close to 90 deg
% 
% V=10;  %the velocity should be given
% 
% mu=V*cos(alpha)/Vt;
% lam=V*sin(alpha)/Vt+a*sig/16*(sqrt(1+24*th0/(a*sig))-1);
% 
% % lam=mu*tan(alpha)+ct/(2*sqrt(mu^2+Lam^2));
% 
% ct=sig*a*0.5*(th0/3*(1+1.5*mu^2)-0.5*lam);
% 
% T=ct*ro*Vt^2*A
