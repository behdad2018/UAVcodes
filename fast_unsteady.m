% simplified unsteady 3D wing lifting line method by B. Davoudi 7/27/2016
tic
clc; clear all;
Nxw=8;Nxt=2;               % chordwise panel per hald wing
Nyw=10;Nyt=Nyw/2;               % spanwise panel per hald wing
ARw=176.5^2/7150;             % Aspect ratio
ARt=3.5;
bw=1;bt=bw*64.4/176.5;         % wing span
%bt=1;ARt=ARw;
trw=1.8;trt=2;              % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoa=degtorad(7);             % angle of attack
uinf=1;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
aoaf_Lw=degtorad(0);aoaf_Lt=degtorad(0);
aoaf_Rw=degtorad(0);aoaf_Rt=degtorad(0);

% geometry
%[x,y,z,xcol,ycol,zcol,n,dl_x,dly,S,alpha,cr]=geometry(AR,b,tr,Nx,Ny,Lam,dih,aoa,aoaf_L,aoaf_R);
[xw,yw,zw,xcolw,ycolw,zcolw,nw,dl_xw,dlyw,Sw,alphaw,crw]=geometry(ARw,bw,trw,Nxw,Nyw,Lam,dih,aoa,aoaf_Lw,aoaf_Rw);
[xt,yt,zt,xcolt,ycolt,zcolt,nt,dl_xt,dlyt,St,alphat,crt]=geometry(ARt,bt,trt,Nxt,Nyt,Lam,dih,aoa,aoaf_Lt,aoaf_Rt);
xt=bw/2+xt;xcolt=bw/2+xcolt;
zt=bw/16+zt;zcolt=bw/16+zcolt;

Nx=Nxw+Nxt;
Ny=Nyw;

x=[xw;[xt,zeros(Nxt,2*Ny-2*Nyt)]];y=[yw;[yt,zeros(Nxt,2*Ny-2*Nyt)]];z=[zw;[zt,zeros(Nxt,2*Ny-2*Nyt)]];
xcol=[xcolw;[xcolt,zeros(Nxt,2*Ny-2*Nyt)]];ycol=[ycolw;[ycolt,zeros(Nxt,2*Ny-2*Nyt)]];
zcol=[zcolw;[zcolt,zeros(Nxt,2*Ny-2*Nyt)]];n=[nw,nt];dl_x=[dl_xw,dl_xt];dly=[dlyw,dlyt];
alpha=[alphaw;[alphat,zeros(Nxt,2*Ny-2*Nyt)]];

%time step
dt=0.5*min(dl_xw)*Nxw/uinf;   %%%%%%%%%%%%%%%%%%%%%% 1/4

%maximum simulation time
tmax=floor(bw/(uinf*dt));

[Gs,Am_wing,A,a]=fast_steady(x,y,z,xcol,zcol,n,dl_x,dly,Nxw,Nxt,Nyw,Nyt,u,alpha,Lam,dih,bw);

%[Gs,Am_wing,a]=fast_steady(x,y,z,xcol,zcol,n,dl_x,dly,Nx,Ny,u,alpha,Lam,dih,b);
% for i=1:tmax
%     xw(i,:)=xt(Nx,:)+dl_x*cos(aoa)+uinf*(i-1)*dt;
%     yw(i,:)=yt(Nx,:);
%     zw(i,:)=zt(Nx,:)-dl_x*sin(aoa)*cos(dih);
% end
for i=1:tmax
    xww(i,:)=xw(Nxw,:)+dl_x(1:2*Nyw).*cos(alphaw(end,:))+uinf*(i-1)*dt;
    yww(i,:)=yw(Nxw,:);
    zww(i,:)=zw(Nxw,:)-dl_x(1:2*Nyw).*sin(alphaw(end,:))*cos(dih);
    if i<=tmax/2   % size of the tail wake
        xwt(i,:)=xt(Nxt,:)+dl_x(2*Nyw+1:2*(Nyt+Nyw)).*cos(alphat(end,:))+uinf*(i-1)*dt;
        ywt(i,:)=yt(Nxt,:);
        zwt(i,:)=zt(Nxt,:)-dl_x(2*Nyw+1:2*(Nyt+Nyw)).*sin(alphat(end,:))*cos(dih);
    end
end
%surf(x,y,z);hold on;
% figure
% surf(xt,yt,zt); hold on;surf(xw,yw,zw);surf(xww,yww,zww);surf(xwt,ywt,zwt);axis('equal');
% figure
% plot(Gs); hold on;
[cLnewvs,cls_section,dps]=force_calc(Nxw,Nyw,a(:,1:2*Nxw*Nyw),Gs(1:2*Nxw*Nyw),Gs(1:2*Nxw*Nyw),1,1,Sw,dl_xw,dlyw,uinf,aoa,Lam,dih,aoaf_Lw,aoaf_Rw);
plot(yw(1,:),cls_section);hold on;xlabel('span');ylabel('C_l');grid minor;
%% unsteady response to a gust

% the initial wake vortex strength obtained by steady solution
Gww=repmat(Gs(1+2*(Nxw-1)*Nyw:2*Nxw*Nyw),tmax,1);Gwt=repmat(Gs(1+2*Nxw*Nyw+2*(Nxt-1)*Nyt:2*Nxt*Nyt+2*Nxw*Nyw),floor(tmax/2),1);
% gust info
vtmax=tmax;
xvor=-bw/2;
yvor=bw/2-bw/32;
zvor=-bw/16;
rc=0.162*bw/8;
% Gv=1.122*uinf*b/8;
Gv=0;

for t=1:vtmax
    t
    xvor=xvor+uinf*dt;
    
    %% total of the panels in the wing and wake
    for i=1:2*Nxw*Nyw+2*Nxt*Nyt
        
        % once the panels on the wing are considered, tail panels are the next
        % starting from 2*Nxw*Nyw+1 to 2*Nxw*Nyw+2*Nxt*Nyt
        
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
        
        %% wake
        %% wing wake
        Nww=tmax;
        wake2vw=0;a2w=0;
        
        for j=1:2*Nww*Nyw
            index1ww=ceil(j/(2*Nyw));
            index2ww=j-(2*Nyw)*(index1ww-1);
            xw1=xww(index1ww,index2ww);      yw1=yww(index1ww,index2ww);      zw1=zww(index1ww,index2ww);
            [a1w,wake1vw]=vortexring(n(:,i),uinf*dt,dly(1),0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,Gww(j));
            wake2vw=wake2vw+wake1vw;
            a2w=a2w+a1w;
        end
        
        %% tail wake
        Nwt=floor(tmax/2);
        wake2vt=0;a2t=0;
        
        for j=1:2*Nwt*Nyt
            index1wt=ceil(j/(2*Nyt));
            index2wt=j-(2*Nyt)*(index1wt-1);
            xw1=xwt(index1wt,index2wt);      yw1=ywt(index1wt,index2wt);      zw1=zwt(index1wt,index2wt);
            [a1t,wake1vt]=vortexring(n(:,i),uinf*dt,dly(2),0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,Gwt(j));
            wake2vt=wake2vt+wake1vt;
            a2t=a2t+a1t;
        end
        
        waken(1,i)=wake2vw+wake2vt;        % total effect of the wakes on the panel i
        %% vortex induced effect
        [Vv,Vvn(i)]=vortexline(n(:,i),xcol1,ycol1,zcol1,xvor,yvor,zvor,xvor,yvor+bw/16,zvor,Gv);
        a(:,i)=a2w+a2t+Vv;
    end
    
    % RHS
    RHSn=dot(-repmat(u',1,2*Nxw*Nyw+2*Nxt*Nyt),n)'-waken'-Vvn';
    % RHSn=dot(-repmat(u',1,2*Nx*Ny),n)'-waken'-reshape(Vv,2*Nx*Ny,1);
    G(:,t)=Am_wing*RHSn;    % vorticity on the wing
    
    % passing the each row of the wake to the next level
    for i=1:tmax-1
        Gww(1+2*i*Nyw:2*(i+1)*Nyw)=Gww(1+2*(i-1)*Nyw:2*i*Nyw) ;
        if i <=tmax/2-1   % size of the tail wake
            Gwt(1+2*i*Nyt:2*(i+1)*Nyt)=Gwt(1+2*(i-1)*Nyt:2*i*Nyt) ;
        end
    end
    
    Gww(1:2*Nyw)= G(1+2*(Nxw-1)*Nyw:2*Nxw*Nyw,t);
    Gwt(1:2*Nyt)= G(1+2*Nxw*Nyw+2*(Nxt-1)*Nyt:2*Nxw*Nyw+2*Nxt*Nyt,t);
    
    [cLnewv(t),cl_section(:,t),dp(:,:,t)]=force_calc(Nxw,Nyw,a(:,1:2*Nxw*Nyw),G(1:2*Nxw*Nyw,:),Gs(1:2*Nxw*Nyw,1),t,dt,Sw,dl_x(1:2*Nyw),dly(1),uinf,aoa,Lam,dih,aoaf_Lw,aoaf_Rw);
    [cLnewvt(t),cl_sectiont(:,t),dpt(:,:,t)]=force_calc(Nxt,Nyt,a(:,1+2*Nxw*Nyw:2*Nxw*Nyw+2*Nxt*Nyt),G(1+2*Nxw*Nyw:2*Nxw*Nyw+2*Nxt*Nyt,:),Gs(1+2*Nxw*Nyw:2*Nxw*Nyw+2*Nxt*Nyt,1),t,dt,St,dl_x(1+2*Nyw:2*(Nyw+Nyt)),dly(2),uinf,aoa,Lam,dih,aoaf_Lt,aoaf_Rt);
end
toc
%% figure
%plot(linspace(0,vtmax*dt*uinf/cr,vtmax),cLnewv);xlabel('time'); ylabel('C_L');grid minor;hold on;
%figure
%plot(y,cl_section(:,10));
%axis([0 9 0 0.55]);
% surf(x,y,z,[G(1:2*Ny,34)';G(2*Ny+1:4*Ny,34)';G(4*Ny+1:6*Ny,34)';G(6*Ny+1:8*Ny,34)'])
figure
surf(repmat(linspace(0,vtmax*dt*uinf/crw,vtmax),2*Ny,1),repmat(y(1,1:2*Ny),vtmax,1)',cl_section)
xlabel('time');ylabel('span');zlabel('c_l');
view(-150,15);

%% plotting real wing -- converting the middle points data to edges.

sp(dp,xw,yw,zw,bw,bt,Nxw,Nyw,trw,Sw,Lam,dih,aoa,aoaf_Lw,aoaf_Rw,tmax,.8,0);
sp(dpt,xt,yt,zt,bt,bw,Nxt,Nyt,trt,St,Lam,dih,aoa,aoaf_Lt,aoaf_Rt,tmax,.8,1);

[X,Y,Z] = cylinder(bw/9);
fill3(repmat(-.2,1,length(X(1,:))),X(1,:),Y(1,:),'g');