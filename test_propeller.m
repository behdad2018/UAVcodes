tic
clc; clear all;
%close all;
%% wing config

Nxw=4; Nxt=0;                 % chordwise panel per hald wing
Nyw=60;Nyt=0/2;               % spanwise panel per hald wing
ARw=176.5^2/7150;             % Aspect ratio
ARt=3.5;
bw=2;
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

%R=0.1*bw;               % propeller radius
%c=0.15*R;                % propeller chord
%rpm=2500;               % rpm

nb=4;                    % number of blade

J=0.7;
R=0.16*bw;               % propeller radius
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
%[xt,yt,zt,xcolt,ycolt,zcolt,nt,dl_xt,dlyt,St,alphat,crt]=geometry(ARt,bt,trt,Nxt,Nyt,Lam,dih,aoat,aoaf_Lt,aoaf_Rt);
%xt=bw/2+xt;xcolt=bw/2+xcolt;
%zt=bw/16+zt;zcolt=bw/16+zcolt;

% Nx=Nxw+Nxt;
% Ny=Nyw;
% 
% x=[xw;[xt,zeros(Nxt,2*Ny-2*Nyt)]];y=[yw;[yt,zeros(Nxt,2*Ny-2*Nyt)]];z=[zw;[zt,zeros(Nxt,2*Ny-2*Nyt)]];
% xcol=[xcolw;[xcolt,zeros(Nxt,2*Ny-2*Nyt)]];ycol=[ycolw;[ycolt,zeros(Nxt,2*Ny-2*Nyt)]];
% zcol=[zcolw;[zcolt,zeros(Nxt,2*Ny-2*Nyt)]];n=[nw,nt];dl_x=[dl_xw,dl_xt];dly=[dlyw,dlyt];
% alpha=[alphaw;[alphat,zeros(Nxt,2*Ny-2*Nyt)]];

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

%%
for i=1:2*Nxw*Nyw+2*Nxt*Nyt
    i;
    % for main wing
    index1col=ceil(i/(2*Nyw));
    index2col=i-(2*Nyw)*(index1col-1);
    
    % once we start consdiering the tail panels
    if index1col>Nxw
        index1col=Nxw+ceil((i-2*Nxw*Nyw)/(2*Nyt));        % goes from Nxw+1 to Nxt
        index2col=i-2*Nxw*Nyw-(2*Nyt)*(index1col-Nxw-1);  % goes from 1 to Nyt
    end
    % collocation points for both wing and tail
    xcol1=xcol(index1col,index2col);
    ycol1=y(index1col,index2col);
    zcol1=zcol(index1col,index2col);
    Vv3=zeros(3,1);Vvn3=0;
    
    %% new way by knowing the induced velocity from the paper
    
    Vvn(i)=0;

    if ycol1<=max(max(max(ys))) && ycol1>=min(min(min(ys)))
        er=abs(abs(ycol1)-abs(yp))-r;
        in=find(abs(er)==min(min(abs(er)))) % this fin the best index for which r(in) and ycol1 are closest
        
        if   r(in)>0.2*r(end)
            if ycol1>yp
                
                
               % Vvn(i)=om*r(in)-nb*gd(in)/(2*pi*r(in))* 0.5;
               Vvn(i)=-nb*gd(in)/(2*pi*r(in));
              %   Vvn(i)=-nb*gd(in)/(4*pi*r(in))*(1+abs(xcol1-xp)/sqrt((xcol1-xp)^2+r(in)^2));
            else
                
              %  Vvn(i)=-om*r(in)+nb*gd(in)/(2*pi*r(in))* 0.5;
               Vvn(i)=nb*gd(in)/(2*pi*r(in));
              %  Vvn(i)=nb*gd(in)/(4*pi*r(in))*(1+abs(xcol1-xp)/sqrt((xcol1-xp)^2+r(in)^2)) ;
            end
            vxx=0.5*om/(2*pi)*nb*gd(in)/uinf*(1+abs(xcol1-xp)/sqrt((xcol1-xp)^2+r(in)^2));
        %    Vx=vxx+uinf;
        end
        % note that we have to subtract uinf from Vx since this only
        % velocity induced by the propeller
       % ap(:,i)=[Vx-uinf;0;Vvn(i)];
       ap(:,i)=ap(:,i)+[vxx*cos(aoaw)*0.5;0;Vvn(i)*cos(aoaw)*0.5-vxx*sin(aoaw)];
        %        ap(:,i)=[Vx-uinf;0;0];
        %        ap(:,i)=[0;0;Vvn(i)];
        % fixing the free stream because the propeller induces velocity
        
      %  urhs(1,i)= Vx;
        urhs(1,i)= uinf;
      %  Rs=R*sqrt((1+vxx/uinf)/(1+vxx*((1+(xcol1-xp)/sqrt((xcol1-xp)^2+R^2)))))
    end
    
    %% old way of doing it
    %     % if abs(ycol1)<max(max(max(ys)))
    %     for ri=1:49      % denotes different layer of the tubes
    %
    %         for k=1:19   % denotes the azimuth angle
    %             dr=r(ri)/50;dth=(2*pi/20);
    %             [Vv1,Vvn1]= vortexline(n(:,i),xcol1,ycol1,zcol1,xp,ys(1,k,ri),zs(1,k,ri),bw+xp,ys(1,k,ri),zs(1,k,ri),gx(ri)*dr*dth);
    %             [Vv2,Vvn2]= vortexline(n(:,i),xcol1,ycol1,zcol1,xp,ys(1,k,ri),zs(1,k,ri),xp,ys(1,k,ri+1),zs(1,k,ri+1),gd(ri)*dth);
    %             %            Vv1=0;Vvn1=0;
    %             %           Vv2=0;Vvn2=0;
    %
    %             for j=1:9 % sweep for stream-wise rings -- arc elements, they have the same x "only" in a ring,
    %
    %                 [Vv3,Vvn3]=vortexline(n(:,i),xcol1,ycol1,zcol1,xs(j,k,ri),ys(j,k,ri),zs(j,k,ri),xs(j,k,ri),ys(j,k+1,ri),zs(j,k+1,ri),gt(ri)*dr/9);
    %                 %v3=0;Vvn3=0;
    %                 Vvn3=Vvn1+Vvn2+Vvn3;
    %                 Vv3=Vv1+Vv2+Vv3;
    %
    %             end
    %         end
    %         % end
    %     end
    %     ap1(:,i)=Vv3;
    %     Vvn(i)=Vvn3;
    
end

%%

RHSn=dot(-urhs,n)'-Vvn'*0.5;
% RHSn=dot(-repmat(u',1,2*Nx*Ny),n)'-waken'-reshape(Vv,2*Nx*Ny,1);
G=A\RHSn;    % vorticity on the wing
RHSn2=dot(-repmat(u',1,2*Nxw*Nyw+2*Nxt*Nyt),n)';
G2=A\RHSn2;    % vorticity on the wing

%%
% figure
% surf(xw,yw,zw);hold on; 
% surf(xt,yt,zt)
% axis('equal');
% surf(xs(:,:,end),ys(:,:,end),zs(:,:,end));

figure
[cLnewvs,cls_section,dps]=force_calc(Nxw,Nyw,ap(:,1:2*Nxw*Nyw),G(1:2*Nxw*Nyw),G(1:2*Nxw*Nyw),1,1,Sw,dl_xw,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw);
[cLnewvs,cls_section2,dps1]=force_calc(Nxw,Nyw,a(:,1:2*Nxw*Nyw),G2(1:2*Nxw*Nyw),G2(1:2*Nxw*Nyw),1,1,Sw,dl_xw,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw);
data=dlmread('-.txt','',[1 0 288 2]);
plot(data(:,1),data(:,3));hold on;plot(data(:,1),data(:,2))
plot(yw(1,:),cls_section);hold on;plot(yw(1,:),cls_section2);xlabel('span');ylabel('C_l');grid minor;legend('prop on','prop off')