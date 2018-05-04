% theta=[-2 0 2 4 8 10 12 14];
% clroot=[0.17 0.36 0.55 0.73 1.07 1.22 1.34 1.41];
% clmid=[0.31 0.506 0.7 0.88 1.1 1.2 1.31 1.128];
% cltip=[0.38 0.31 0.58 0.81 .88 0.92 .95 1.0];
%
% plot(theta,clroot);hold on
% plot(theta,clmid);plot(theta,cltip);grid minor;legend('root','mid','tip');
%
% hold on;grid minor;
% exp=load('fort.1100_3d');
% plot(exp(:,1),exp(:,2));cm
% ylim([0 1.5]);

% to generate axial and swirl plot for aiaa paper 2017

clear all, close all; 
clc;
bw=2;
uinf=50 ;
nb=4;                    % number of blade
J=0.7;
R=0.16*bw;               % propeller radius
Rs=R;
ARp=7;
Ap=Rs^2/ARp;
c1=[0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.989, 0.994, 0.988, 0.961, 0.894, 0.732, 0.513];

%cm=Ap/trapz(linspace(0,1,20),c1);
cm=Ap/trapz(linspace(0,1,20)*Rs,c1)

rpm=uinf/J/(2*R)*60 ;        % rpm
ro=1;                   % density
uinfp=uinf;                 % free stream velocity

A=pi*R^2;               % disk area
sig=nb*cm/(pi*R);        % solidity
om=rpm*(2*pi/60);       % rad/second
gam=ro*2*pi*cm*R^4/1;    % lock number
Vt=om*R;                % tip velocity
alpha_p=degtorad(90);   % for a propeller it is close to 90 deg

%xp=-0.20;yp=0;zp=0;
xp=0;yp=0;zp=0;
num1=8;                 % this is the number of sections in each circle
num2=16;                 % number of rings (sqrt must be an integer)
rnum=20;                 % number of cylinders

[xs,ys,zs,r1,gx,gt,gd,cl_asitav]=propeller(bw,xp,yp,zp,num1,num2,rnum,nb,rpm,R,cm,uinfp,alpha_p,uinfp);
r=r1*R;
location=[0.14 0.47 1.72 3]*R;
 c_a=0.65;
%c_a=1;
c_s=.9;

Vvn=zeros(20,4);vxx=zeros(20,4);
for i=1:4
    
     
    vxa=c_a*  0.5*om/(2*pi)*nb*gd/uinf;     % axial velocity increase at the disk
    vxaa=mean(vxa);
    % finding slip-stream radius in the downstream
    
    Rs(i)=R*sqrt((1+vxaa/uinf)/(1+vxaa/uinf*((1+abs(location(i)-xp)/sqrt((location(i)-xp)^2+R^2)))))
    uinf
    
    Vvn(:,i)=c_s*nb*gd./(2*pi*r*Rs(i)/R);
    vxx=vxa.*(1+abs(location(i)-xp)./sqrt((location(i)-xp)^2+(r*Rs(i)/R).^2));
    Vx(:,i)=uinf+vxx;  
        
    rs(:,i)=r(4:end)*Rs(i)/R^2;
    swirl(:,i+1)=atan(Vvn(4:end,i)./Vx(4:end,i))*180/pi;
     
    
    figure(1)
    plot(rs(:,i),Vx(4:end,i)/uinf);hold on;legend('x= 0.14R','x= 0.47R','x=1.72','x=3R');
     
     figure(2)
     plot(rs(:,i),swirl(:,i+1));hold on;xlabel('r/R');ylabel('swirl angle (deg)');hold on
     legend('x= 0.14R','x= 0.47R','x=1.72','x=3R')
end

%  save('Jaiaa_axial.txt','text','-ascii')