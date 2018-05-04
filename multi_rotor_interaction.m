% multi-rotor prop-wing interaction
clear all; clc;

%% wing config

Nxw=4; Nxt=0;                 % chordwise panel per hald wing
Nyw=61;Nyt=0/2;               % spanwise panel per hald wing
ARw=4;             % Aspect ratio
ARt=3.5;
bw=8;
bt=bw*64.4/176.5;          % wing span
bt=1;ARt=ARw;
trw=1;trt=2;              % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoaw=degtorad(5);             % angle of attack of the wing
aoat=degtorad(0);             % angle of attack of the tale
uinf=50 ;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
urhs=repmat(u',1,2*Nxw*Nyw+2*Nxt*Nyt);  % this will be modified to the effect of the propeller
aoaf_Lw=degtorad(0);aoaf_Lt=degtorad(0);
aoaf_Rw=degtorad(0);aoaf_Rt=degtorad(0);

%% propeller config

nb=4;                    % number of blade
J=0.7;
R=0.04*bw;               % propeller radius
Rs=R;
ARp=7;
Ap=(0.32*bw)^2/ARp;
c1=[0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.989, 0.994, 0.988, 0.961, 0.894, 0.732, 0.513];
cm=Ap/trapz(linspace(0,1,20),c1);

rpm=uinf/J/(2*R)*60 ;        % rpm
ro=1;                   % density
uinfp=uinf;                 % free stream velocity

ARw=4;             % Aspect ratio

A=pi*R^2;               % disk area
sig=nb*cm/(pi*R);        % solidity
om=rpm*(2*pi/60);       % rad/second
gam=ro*2*pi*cm*R^4/1;    % lock number
Vt=om*R;                % tip velocity
alpha_p=degtorad(90);   % for a propeller it is close to 90 deg


yprop=[-3,-1.5,0,1.5,3];

%% location of the propeller

%xp=-0.20;yp=0;zp=0;
xp=-0.25-0.435;yp=0;zp=0;
num1=8;                 % this is the number of sections in each circle
num2=16;                 % number of rings (sqrt must be an integer)
rnum=20;                 % number of cylinders
[xs,ys,zs,r1,gx,gt,gd,cl_asitav]=propeller(bw,xp,yp,zp,num1,num2,rnum,nb,rpm,R,cm,uinfp,alpha_p,uinfp);
r=r1*R;

%% lifting surfaces locations and geometry

[xw,yw,zw,xcolw,ycolw,zcolw,nw,dl_xw,dlyw,Sw,alphaw,crw]=geometry(ARw,bw,trw,Nxw,Nyw,Lam,dih,aoaw,aoaf_Lw,aoaf_Rw);

Nx=Nxw+Nxt;
Ny=Nyw;
%
x=xw;y=yw;z=zw;
xcol=xcolw;ycol=ycolw;
zcol=zcolw;n=nw;dl_x=dl_xw;dly=dlyw;
alpha=alphaw;


[Gs,Am_wing,A,a]=fast_steady(x,y,z,xcol,zcol,n,dl_x,dly,Nxw,Nxt,Nyw,Nyt,u,alpha,Lam,dih,bw);
% propeller induced vel will be over written on this
ap=a;

% coefficients to fix axial and swirl velocity obtained from prop alone study
c_a=0.7 * .5;
c_s=0.9 * .5;
%%
for i=1:2*Nxw*Nyw
    
    Vvn(i)=0;
    
    % for main wing
    index1col=ceil(i/(2*Nyw));
    index2col=i-(2*Nyw)*(index1col-1);
    
    % collocation points for both wing and tail
    xcol1=xcol(index1col,index2col);
    ycol1=y(index1col,index2col);
    zcol1=zcol(index1col,index2col);
    
    % finding slip-stream radius in the downstream
    
    vxa=c_a*  0.5*om/(2*pi)*nb*gd/uinf;     % axial velocity increase at the disk
    vxaa=mean(vxa);
    % Rs(index1col)=R*sqrt((1+vxaa/uinf)/(1+vxaa/uinf*((1+abs(xcol1-xp)/sqrt((xcol1-xp)^2+R^2)))));
    Rs(index1col)=R;
    p=find(abs(yprop-ycol1)==min(abs(yprop-ycol1)))
    
    if ycol1<=Rs(index1col)+yprop(p) && ycol1>=-Rs(index1col)+yprop(p)
        er=abs(abs(ycol1)-abs(yprop(p)))-r*Rs(index1col)/R;
        in=find(abs(er)==min(min(abs(er)))); % this fin the best index for which r(in) and ycol1 are closest
        
        if   r(in)>0.2*R
            
            if ycol1>yprop(p)
                Vvn(i)=c_s* -nb*gd(in)/(2*pi*r(in)*Rs(index1col)/R);
            else
                Vvn(i)=c_s* nb*gd(in)/(2*pi*r(in)*Rs(index1col)/R);
            end
            vxx=vxa(in)*(1+abs(xcol1-xp)/sqrt((xcol1-xp)^2+(r(in)*Rs(index1col)/R)^2));
            ap(:,i)=ap(:,i)+[vxx*cos(aoaw);0;Vvn(i)*cos(aoaw)-vxx*sin(aoaw)];
%         else
%              urhs(:,i)=0;
%              % where the naccele is the lift is killed
%              ap(:,i)=[-50;0;0];
%              a(:,i)=[-50;0;0];
        end
        if or(and(r(in)<=0.35*R,index1col==1),and(r(in)<0.3*R,index1col==2))        
            urhs(:,i)=0;
            % where the naccele is the lift is killed
            ap(:,i)=[-50;0;0];
            a(:,i)=[-50;0;0];
        elseif or(and(r(in)<=0.25*R,index1col==3),and(r(in)<=0.2*R,index1col==4))
            urhs(:,i)=0;
            % where the naccele is the lift is killed
            ap(:,i)=[-50;0;0];
            a(:,i)=[-50;0;0];
        end        
    end
end

%%

RHSn=dot(-urhs,n)'-Vvn';
G=A\RHSn;    % vorticity on the wing
RHSn2=dot(-urhs,n)';
G2=A\RHSn2;    % vorticity on the wing

%%
figure
%% 



[cLnewvs,cls_section,dps]=force_calc(Nxw,Nyw,ap(:,1:2*Nxw*Nyw),G(1:2*Nxw*Nyw),G(1:2*Nxw*Nyw),1,1,Sw,dl_xw,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw);
[cLnewvs,cls_section2,dps1]=force_calc(Nxw,Nyw,a(:,1:2*Nxw*Nyw),G2(1:2*Nxw*Nyw),G2(1:2*Nxw*Nyw),1,1,Sw,dl_xw,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw);
data=dlmread('-.txt','',[1 0 288 2]);

cls_section(cls_section == 0) = NaN;
cls_section2(cls_section2 == 0) = NaN;


% plot(data(:,1),data(:,3));hold on;plot(data(:,1),data(:,2))
plot(yw(1,:),cls_section);hold on;plot(yw(1,:),cls_section2);xlabel('span');ylabel('C_l');grid minor;
legend('prop on, CFD','prop off, CFD','prop on, model','prop off, model');