function [T,Q,P]=engine(rpm,alpha,beta,gamma,Vx,Vy,Vz)

% Individual engine code a typical quad-copter (only one of the prop). This
% code provides Thrust, Torque and Power generated/required for a given
% flight condition variables and copter geometric characterstics -- written by B.
% Davoudi, 2017 - Computational Aerosciences Laboratory at University of
% Michigan

tic
nb=2;                     % number of blade
rnum =20;                 % number of radial points on a blade
%R=0.12;                    % propeller radius
R=3*0.0254;
cm=0.08;                  % maxium chord  
rho=1.15;                 % air density
% rpm=5000;               % rpm
% R=0.2;                  % propeller radius [m]
% c=0.05;                 % propeller chord [m]
% Vx,Vy,Vz                % the relative velocity of the copter wrt wind
% alpha, beta and gamma   % the angles of prop axis wrt ground x,y,and axes in rad

A=pi*R^2;                % disk area
sig=nb*cm/(pi*R);        % solidity
om=rpm*2*pi/60;          % rad per second
Vt=om*R;                 % tip velocity

prop_axis=[cos(alpha),cos(beta),cos(gamma)];
rel_wind=[Vx,Vy,Vz];
V=dot(prop_axis,rel_wind);

lamc=V*cos(gamma)/Vt                        % climb inflow ratio
r=linspace(1/rnum,1,rnum);                    % propeller radial

th0=[6.3,1];
r0=[0,1];the0=interp1(r0,th0,r);
cl0=[1.5*pi 1.7*pi];
cla=interp1(r0,cl0,r);

th=degtorad(the0+[56.7,56.7,56.7,56.7, 51.3, 47.0, 43.5, 40.7, 38.3, 36.3, 34.5, 32.8, 31.3, 29.8, 28.3, 26.9, 25.5, 24.0, 22.4, 20.8]);
th=degtorad(the0+[56.7,56.7,56.7,56.7, 51.3, 47.0, 43.5, 40.7, 38.3, 36.3, 34.5, 32.8, 31.3, 29.8, 28.3, 26.9, 25.5, 24.0, 22.4, 20.8]);
c=cm*[0, 0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.989, 0.994, 0.988, 0.961, 0.894, 0.732, 0.513];

%% BEMT to be used to find the  / bound vortex on the disk

% Prandtl's tip-loss function
r(end)=0.99;r(1)=0.01;
lam=zeros(1,rnum);
phi=lam;
F=phi;
for i=1:length(r)
    Fold=1;dF=1;
    while abs(dF)>.1
        
        lam(i)=sqrt((sig*cla(i)/(16*Fold)-lamc*0.5)^2 + sig.*cla(i)*th(i)*r(i)*0.125/Fold)- (sig*cla(i)/(16*Fold)-lamc*0.5);
        phi(i)=atan(lam(i)/r(i));
        froot=0.5*nb*r(i)/((1-r(i))*phi(i));
        ftip=0.5*nb*(1-r(i))/(r(i)*phi(i));
        
        %% based on pp 142 and 144 of Leishman -- does not work
        %  Fnew=(2/pi)*acos(exp(-froot*ftip))
        
        %% My propossal  
        Fnew=(2/pi)*acos(exp(-froot))*(2/pi)*acos(exp(-ftip));
        
        %%
        dF=abs(Fnew-Fold);
        Fold=Fnew;         
    end
   
    F(i)=Fnew;
end

cl=cla.*(th-phi);

% thrust calculations

L=cl.*0.5*1.225.*sqrt((r*R*om).^2+V^2).^2.*c;
T=trapz(r*R,L)*nb;
% plot(r*R,L);

% exterted forced in x,y, and directions
% Tx=T*cos(alpha)
% Ty=T*cos(beta)
% Tz=T*cos(gamma)

%% Torque  -- July 2017

% blade profile drag
% Cd0=0.0081-0.0216*aoa+0.4*aoa^2; from helicopter literature, note aoa is the angle of attack - page 14, chapter 2 of Friedmann notes
Cd0=0.012;

% thrust coefficient
Ct=T/(rho * pi * R^2 * Vt^2);

% f is the equivalent flat plate area that copter frame occupies in space
% -- note that a quarter of the area should be used since this is model for
% one of the prop, for instance that area is 1 m^2

f= 0.25 *.5;
f=0.25 * 0.3*0.3;
% mu is the advance ratio

mu= sqrt(Vx^2+Vy^2) / Vt

% total power required based on page 80, chapter 2 of Friedmann notes

Cp = sig*Cd0*0.125*(1+4.6*mu^2) + Ct^2/(2*mu) + f * mu^3 /(2*pi*R^2) + Ct*lamc 

Cq = Cp;

Q=Cq * rho * pi * R^2 * Vt^3 / om;
P=Q*om

toc
