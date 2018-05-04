% simplified unsteady 3D wing lifting line method by B. Davoudi 7/27/2016
% used for aiaa paper 2017
tic
clc; clear all;
Nxw=4;                       % chordwise panel per hald wing
Nyw=85;                      % spanwise panel per hald wing
AR=10.75;                    % Aspect ratio
bw=64.5*2*0.0254;            % wing span
tr=1;                        % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoaw=degtorad(7);             % angle of attack of the wing
aoat=degtorad(0);             % angle of attack of the tale
uinf=85.34;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
aoaf_Lw=degtorad(0);
aoaf_Rw=degtorad(0);
urhs=repmat(u',1,2*Nxw*Nyw);

cfd=load('fort.1101');
exper=[8	0.651008
30	0.700044
42	0.656122
43	0.655377
44	0.647852
45	0.641056
46.5	0.631097
48	0.615122
49	0.604804
50	0.595209
51	0.581881
52.5	0.568698
54	0.543273
55	0.527745
56	0.506622
57	0.486198
58.5	0.450677
60.5	0.386243
61.75	0.332016
63	0.266755
64	0.285729];


% geometry
[x,y,z,xcol,ycol,zcol,n,dl_xw,dlyw,Sw,alphaw,cr]=geometry(AR,bw,tr,Nxw,Nyw,Lam,dih,aoaw,aoaf_Lw,aoaf_Rw);

%[x1,y1,z1,xcol1,ycol1,zcol1,n1,dl_x1,dly1,S1,alpha1,cr1]=geometry(5.166,46.5*2,tr,Nxw,Nyw,Lam,dih,4,aoaf_L,aoaf_R);

% alpha(:,20)=0;alpha(:,21)=0;

[Gs,Am_wing,A,a]=fast_steady(x,y,z,xcol,zcol,n,dl_xw,dlyw,Nxw,0,Nyw,0,u,alphaw,Lam,dih,bw);
%plot(G15); hold on;

%[cLnewvs,cls_section,dps]=force_calc(Nxw,Nyw,a,Gs,Gs,1,1,S,dl_x,dly,uinf,aoaw,Lam,dih,aoaf_L,aoaf_R);
%% VG
p=0;

for delz=[0 -2]
% delz=0; % direct hit
cvg=18*0.0254;bvg=(52.5+delz)*0.0254;
ARvg=bvg/cvg;
[xvg,yvg,zvg,xcolvg,ycolvg,zcolvg,nvg,dl_xwvg,dlywvg,Swvg,alphawvg,crvg]=geometry(ARvg,bvg,tr,Nxw,Nyw,0,0,degtorad(8),0,0);
[Gsvg,Am_wingvg,Avg,avg]=fast_steady(xvg,yvg,zvg,xcolvg,zcolvg,nvg,dl_xwvg,dlywvg,Nxw,0,Nyw,0,u,alphawvg,0,0,bvg);

%%
% plot(exper(:,1),exper(:,2),'+');hold on;grid on;
% plot(cfd(:,1)*64.5,cfd(:,2)*1.065)
% plot(y(1,floor(end/2)+1:end)/0.0254,cls_section(floor(end/2)+1:end)*1.065);
% legend('Experiment','CFD','Vortex Method');

%%
% used for different span-wise location at delz=0
%for yv=[46.5 52.5 58.5]*0.0254
 p=p+1;  
for i=1:2*Nxw*Nyw
    Vvn(i)=0;
    % for main wing
    index1col=ceil(i/(2*Nyw));
    index2col=i-(2*Nyw)*(index1col-1);

    % collocation points for both wing and tail
    xcol1=xcol(index1col,index2col);
    ycol1=y(index1col,index2col);
    zcol1=zcol(index1col,index2col);
    
    %% induced velocity from the vortex
    yv=46.5*0.0254;
    % bw/2*0.75;     % location of spot the vortex hit
    zv=delz*0.0254;
    xv=1;
    %rc=0.02;
    rc=1.5*0.0254;  % based on the experimentalist recommendation
    r=sqrt((yv-ycol1)^2+zv^2);        % propeller radial   
    gam=3.5;
    gam=mean([Gsvg(Nyw) Gsvg(3*Nyw) Gsvg(5*Nyw) Gsvg(7*Nyw)]);
   if ycol1>yv
    Vvn(i)=gam/(2*pi*r)*(1-exp(-r^2/rc^2));
   else
    Vvn(i)=-gam/(2*pi*r)*(1-exp(-r^2/rc^2));
   end
    ap(:,i)=[0;0;Vvn(i)*cos(aoaw)];
    urhs(1,i)= uinf;
end

%%

RHSn=dot(-urhs,n)'-Vvn';
% RHSn=dot(-repmat(u',1,2*Nx*Ny),n)'-waken'-reshape(Vv,2*Nx*Ny,1);
G=A\RHSn;    % vorticity on the wing
RHSn2=dot(-repmat(u',1,2*Nxw*Nyw),n)';
G2=A\RHSn2;    % vorticity on the wing
%figure
[cLnewvs,cls_section,dps]=force_calc(Nxw,Nyw,ap(:,1:2*Nxw*Nyw),G(1:2*Nxw*Nyw),G(1:2*Nxw*Nyw),1,1,Sw,dl_xw,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw);
[cLnewvs,cls_section2,dps1]=force_calc(Nxw,Nyw,a(:,1:2*Nxw*Nyw),G2(1:2*Nxw*Nyw),G2(1:2*Nxw*Nyw),1,1,Sw,dl_xw,dlyw,uinf,aoaw,Lam,dih,aoaf_Lw,aoaf_Rw);
%plot(y(1,:)/0.0254,cls_section*1.065);hold on;plot(y(1,:)/0.0254,cls_section2*1.065);xlabel('span');ylabel('C_l');grid minor;

aiaa_cls_section(:,1)=y(1,end/2+1:end)/0.0254;
aiaa_cls_section(:,p+1)=cls_section(end/2+1:end);

% aiaa_cls_section(:,5)=cls_section2(end/2+1:end);

plot(aiaa_cls_section(:,1),aiaa_cls_section(:,p+1));hold on;xlabel('span');ylabel('C_l');grid minor;
%plot(aiaa_cls_section(:,1),aiaa_cls_section(:,5));
end
