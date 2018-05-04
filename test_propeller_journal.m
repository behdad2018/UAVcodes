tic
clc; clear all; close all;
%close all;
%% wing config

Nxw=15; Nxt=0;                 % chordwise panel per hald wing
Nyw=30;Nyt=0/2;               % spanwise panel per hald wing
ARw=4;             % Aspect ratio
ARt=3.5;
bw=2;
bt=bw*64.4/176.5;          % wing span
bt=1;ARt=ARw;
trw=1;trt=1;              % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoaw=degtorad(5);             % angle of attack of the wing
aoat=degtorad(0);             % angle of attack of the tale
uinf=50 ;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
aoaf_Lw=degtorad(0);aoaf_Lt=degtorad(0);
aoaf_Rw=degtorad(0);aoaf_Rt=degtorad(0);

%% propeller config
%
% nb=4;                    % number of blade
% J=0.7;
% R=0.16*bw;               % propeller radius
% Rs=R;
% ARp=7;
% Ap=(0.32*bw)^2/ARp;

bw=2;                    % wing span of th wing
uinf=50 ;                % cruise velocity
nb=4;                    % number of blade
J=0.7;                   %
R=0.16*bw;               % propeller radius
Rs=R;
ARb=7;                   % Aspect raio of the blade
Ap=Rs^2/ARb;

c1=[0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.989, 0.994, 0.988, 0.961, 0.894, 0.732, 0.513];
cm=Ap/trapz(linspace(0,1,20)*Rs,c1);



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

Nx=Nxw+Nxt;
Ny=Nyw;
%
x=xw;y=yw;z=zw;
xcol=xcolw;ycol=ycolw;
zcol=zcolw;n=nw;dl_x=dl_xw;dly=dlyw;
alpha=alphaw;

%% removing the middle part of the wing due to nacelle
% the technique used here is really elegant but not a good physical model
% bc there is no flow through condition bu in this method flow can pass
% through between the wing

% coefficients to fix axial and swirl velocity obtained from prop alone study
c_a=0.65 * .5;
c_s=0.9 * .5;


% m=0;
% for i=1:2*Nyw 
%     
%     % for main wing
%     index1col=ceil(i/(2*Nyw));
%     index2col=i-(2*Nyw)*(index1col-1);
%     
%     % collocation points for both wing and tail
%     xcol1=xcol(index1col,index2col);
%     ycol1=y(index1col,index2col);
%     zcol1=zcol(index1col,index2col);
%     
%     % finding slip-stream radius in the downstream
%     
%     vxa=c_a*  0.5*om/(2*pi)*nb*gd/uinf;     % axial velocity increase at the disk
%     vxaa=mean(vxa);
%     
%     Rs(index1col)=R*sqrt((1+vxaa/uinf)/(1+vxaa/uinf*((1+abs(xcol1-xp)/sqrt((xcol1-xp)^2+R^2)))));
%     
%     % Rs(index1col)=R;
%     if ycol1<=Rs(index1col) && ycol1>=-Rs(index1col)
%         er=abs(abs(ycol1)-abs(yp))-r*Rs(index1col)/R;
%         in=find(abs(er)==min(min(abs(er)))); % this fin the best index for which r(in) and ycol1 are closest
%         
%         if r(in)<=0.35*R
%             m=m+1;
%             index(m)=index2col;
%         end   
%     end
% end
% 
% 
% x(:,index)=[];
% y(:,index)=[];
% z(:,index)=[];
% 
% xcol(:,index)=[];
% ycol(:,index)=[];
% zcol(:,index)=[];
% alpha(:,index)=[];
% dl_x(index)=[];
% Nyw=Ny-length(index)/2;
% 
% for i=1:Nxw
% remove(1+(length(index)*(i-1)):i*length(index))=2*(i-1)*Ny+index;
% end 
% 
% n(:,remove)=[];

%% steady solution             

[Gs,Am_wing,A,a,~,~]=fast_steady(x,y,z,xcol,zcol,n,dl_x,dly,Nxw,Nxt,Nyw,Nyt,u,alpha,Lam,dih,bw);
% propeller induced vel will be over written on this
ap=a;

% coefficients to fix axial and swirl velocity obtained from prop alone study
% c_a=0.65 * .5;
% c_s=0.9 * .5;
%%
m=0;
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
    Rs(index1col)=R*sqrt((1+vxaa/uinf)/(1+vxaa/uinf*((1+abs(xcol1-xp)/sqrt((xcol1-xp)^2+R^2)))));
    % Rs(index1col)=R;
    if ycol1<=Rs(index1col) && ycol1>=-Rs(index1col)
        er=abs(abs(ycol1)-abs(yp))-r*Rs(index1col)/R;
        in=find(abs(er)==min(min(abs(er)))); % this fin the best index for which r(in) and ycol1 are closest
        
        % spots
        
        if   r(in)>0.3*Rs
            if ycol1>yp
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
        else
            % no induced velocity where nacelle exists
            Vvn(i)=0;
        end
        % Thom Nacelle
        %         if or(and(r(in)<=0.35*R,0<=round(index1col/Nx*100)<25),and(r(in)<0.3*R,25<=round(index1col/Nx*100)<50))
        %             urhs(:,i)=0;
        %             % where the naccele is the lift is killed
        %             ap(:,i)=[-50;0;0];
        %             a(:,i)=[-50;0;0];
        %         elseif or(and(r(in)<=0.25*R,50<=round(index1col/Nx*100)<=75),and(r(in)<=0.2*R,75<=round(index1col/Nx*100)<=100))
        %             urhs(:,i)=0;
        %             % where the naccele is the lift is killed
        %             ap(:,i)=[-50;0;0];
        %             a(:,i)=[-50;0;0];
        %         end
        
        % (no wing here!)
        % reducing the size of matrices and vectors
        if r(in)<=0.3*R
            m=m+1;
            index(m)=i;
        end
    end
end


% reducing the size of matrices and vectors

% A(index,:)=[];
% A(:,index)=[];
% urhs(:,index)=[];
% Vvn(index)=[];
% n(:, index)=[];
% ap(:,index)=[];
% a(:,index)=[];

%%
urhs=repmat(u',1,2*Nx*Nyw);
RHSn=dot(-urhs,n)'-Vvn';
G=A\RHSn;    % vorticity on the wing
RHSn2=dot(-urhs,n)';
G2=A\RHSn2;    % vorticity on the wing

% reconstruct the arrays by assigning zero values at the point there was
% no lifting surface
% 
% m=1;
% for i=1:2*Nx*Ny
%     
%     G3(i,1)=G(i-(m-1));
%     G4(i,1)=G2(i-(m-1));
%     ap3(:,i)=ap(:,i-(m-1));
%     a4(:,i)=a(:,i-(m-1));
%     
%     if m<=length(index)
%         if  i==index(m)
%             G3(i,1)=0;
%             G4(i,1)=0;
%             ap3(:,i)=[0;0;0];
%             a4(:,i)=[0;0;0];           
%             m=m+1
%         end
%     end
%     
% end


%%
figure
%%
[cLnewvs,cls_section,~,~,~,~,~]=force_calc(Nxw,Nyw,ap(:,:),ap(:,:),ap(:,:),bw,G(1:2*Nxw*Nyw),G(1:2*Nxw*Nyw),1,1,Sw,dl_x,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw)
%                                        (Nx,Ny,   a,              a_d,              w_ind_drag,       b, G,             G5,            t,dt,S,dl_x,dly,  uinf,aoa,Lam,dih,aoaf_L,    aoaf_R)
[cLnewvs,cls_section2,~,~,~,~,~]=force_calc(Nxw,Nyw,a(:,1:2*Nxw*Nyw),a(:,1:2*Nxw*Nyw),a(:,1:2*Nxw*Nyw),bw,G2(1:2*Nxw*Nyw),G2(1:2*Nxw*Nyw),1,1,Sw,dl_x,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw);
data=dlmread('-.txt','',[1 0 288 2]);
plot(data(:,1),data(:,3));hold on;plot(data(:,1),data(:,2))
plot(y(1,:),cls_section);hold on;plot(y(1,:),cls_section2);xlabel('span');ylabel('C_l');grid minor;
legend('prop on, CFD','prop off, CFD','prop on, model','prop off, model');



savedata=[y(1,:)',cls_section',cls_section2'];
          

save('Jaiaa_lift_dist_review.txt','savedata','-ascii')


