% to generate axial and swirl plot for aiaa paper 2017
clear all;
close all; 
clc;

bw=2;                    % wing span of th wing
uinf=50 ;                % cruise velocity
nb=4;                    % number of blade
J=0.7;                   % 
R=0.16*bw;               % propeller radius
Rs=R;
ARb=7;                   % Aspect raio of the blade
Ab=Rs^2/ARb;

c1=[0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.989, 0.994, 0.988, 0.961, 0.894, 0.732, 0.513];
cm=Ab/trapz(linspace(0,1,20)*Rs,c1);


%%
rpm=uinf/J/(2*R)*60 ;        % rpm
ro=1;                        % density

A=pi*R^2;                % disk area
sig=nb*cm/(pi*R);        % solidity
om=rpm*(2*pi/60);        % rad/second
gam=ro*2*pi*cm*R^4/1;    % lock number
Vt=om*R;                 % tip velocity
alpha_p=degtorad(90);    % for a propeller it is close to 90 deg


rnum=20;                 % number of radia locations

%%%%%%%%%%%%%%%%%%%%% propeller code %%%%%%%%%%%%%%%%%%%%%%%%%%

% gd is gamma distribution

[r1,G,cl_asitav]=propeller_patricia(rnum,nb,cm,rpm,R,alpha_p,uinf);
                                    
r=r1*R;

% lcations behind the prpeller
location=[0.14 0.47 1.72 3]*R;
xp=0;yp=0;zp=0;

c_a=0.65;
c_s=.9;

Vvn=zeros(20,4);vxx=zeros(20,4);

for i=1:4
    
    vxa=c_a*  0.5*om/(2*pi)*nb*G/uinf;     % axial velocity increase at the disk
    vxaa=mean(vxa);
    
    % finding slip-stream radius in the downstream -- CONTRACTION
    
    Rs(i)=R*sqrt((1+vxaa/uinf)/(1+vxaa/uinf*((1+abs(location(i)-xp)/sqrt((location(i)-xp)^2+R^2)))))
 
    % normal velocity 
    Vvn(:,i)=c_s*nb*G./(2*pi*r*Rs(i)/R);
    
    
    vxx=vxa.*(1+abs(location(i)-xp)./sqrt((location(i)-xp)^2+(r*Rs(i)/R).^2));
    Vx(:,i)=uinf+vxx;  
        
    rs(:,i)=r(4:end)*Rs(i)/R^2;
    
    swirl(:,i)=atan(Vvn(4:end,i)./Vx(4:end,i))*180/pi;

    figure(1)
    plot(rs(:,i),Vx(4:end,i)/uinf);hold on;legend('x= 0.14R','x= 0.47R','x=1.72','x=3R');
 
     figure(2)
     plot(rs(:,i),Vvn(4:end,i)/uinf);hold on;xlabel('r/R');ylabel('rnormal velocity');hold on
     legend('x= 0.14R','x= 0.47R','x=1.72','x=3R')
    
     figure(3)
     plot(rs(:,i),swirl(:,i));hold on;xlabel('r/R');ylabel('swirl angle (deg)');hold on
     legend('x= 0.14R','x= 0.47R','x=1.72','x=3R')
 
    
     
end

vx_save(:,:)=[rs(:,1),Vx(4:end,1)/uinf,rs(:,2),Vx(4:end,2)/uinf...
              rs(:,3),Vx(4:end,3)/uinf,rs(:,4),Vx(4:end,4)/uinf];
          
swirl_save(:,:)=[rs(:,1),swirl(:,1),rs(:,2),swirl(:,2)...
              rs(:,3),swirl(:,3),rs(:,4),swirl(:,4)];          

figure(4)
plot(r,G),ylabel('circulation on the balde');xlabel('r')

save('Jaiaa_axial_review.txt','vx_save','-ascii')
save('Jaiaa_swirl_angle_review','swirl_save','-ascii')
%  save('Jaiaa_axial.txt','text','-ascii')