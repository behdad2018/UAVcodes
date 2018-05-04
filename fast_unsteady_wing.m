% simplified unsteady 3D wing lifting line method by B. Davoudi 7/27/2016
tic
clc; clear all;
Nx=4               % chordwise panel per hald wing
Ny=6               % spanwise panel per hald wing
AR=8;              % Aspect ratio
b=1;               % wing span 
tr=1.8             % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoa=degtorad(5);             % angle of attack
uinf=1;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
aoaf_L=degtorad(0);
aoaf_R=degtorad(0);

% geometry
[x,y,z,xcol,ycol,zcol,n,dl_x,dly,S,alpha,cr]=geometry(AR,b,tr,Nx,Ny,Lam,dih,aoa,aoaf_L,aoaf_R);

%time step
dt=min(dl_x)*Nx/uinf;   %%%%%%%%%%%%%%%%%%%%%% 1/4

%maximum simulation time
tmax=floor(b/(uinf*dt));

[Gs,Am_wing,a]=fast_steady_wing(x,y,z,xcol,zcol,n,dl_x,dly,Nx,Ny,u,alpha,Lam,dih,b);

for i=1:tmax
     xw(i,:)=x(Nx,:)+dl_x*cos(aoa)+uinf*(i-1)*dt;
     yw(i,:)=y(Nx,:);
     zw(i,:)=z(Nx,:)-dl_x*sin(aoa)*cos(dih);
end

%surf(x,y,z);hold on;surf(xt,yt,zt); surf(xw,yw,zw);surf(xwt,ywt,zwt);axis('equal');
% plot(Gs); hold on;
% [cLnewvs,cls_section,dps]=force_calc(Nx,Ny,a,Gs,Gs,1,1,S,dl_x,dly,uinf,aoa,Lam,dih,aoaf_L,aoaf_R);

%%

Gw=repmat(Gs(end-2*Ny+1:end),tmax,1);
vtmax=tmax;   %%%%%%%%%%
xvor=-b/2;

for t=1:vtmax
    t
    xvor=xvor+uinf*dt;
    yvor=b/2-b/32;
    zvor=-b/16;
    rc=0.162*b/8;
    % Gv=1.122*uinf*b/8;
    Gv=0;
    for i=1:2*Nx*Ny
        index1col=ceil(i/(2*Ny));
        index2col=i-(2*Ny)*(index1col-1);
        xcol1=xcol(index1col,index2col);
        ycol1=y(index1col,index2col);
        zcol1=zcol(index1col,index2col);
        % panels on the wake on the pannel on the wing
        Nw=tmax;
        wake2v=0;a2=0;
        
        for j=1:2*Nw*Ny
            index1w=ceil(j/(2*Ny));
            index2w=j-(2*Ny)*(index1w-1);
            xw1=xw(index1w,index2w);      yw1=yw(index1w,index2w);      zw1=zw(index1w,index2w);
            [a1,wake1v]=vortexring(n(:,i),uinf*dt,dly,0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,Gw(j));
            wake2v=wake2v+wake1v;
            a2=a2+a1;
        end
        waken(1,i)=wake2v;
        [Vv,Vvn(i)]=vortexline(n(:,i),xcol1,ycol1,zcol1,xvor,yvor,zvor,xvor,yvor+b/16,zvor,Gv);
        a(:,i)=a2+Vv;
    end
    % RHS
    RHSn=dot(-repmat(u',1,2*Nx*Ny),n)'-waken'-Vvn';
    % RHSn=dot(-repmat(u',1,2*Nx*Ny),n)'-waken'-reshape(Vv,2*Nx*Ny,1);
    G(:,t)=Am_wing*RHSn;    % vorticity on the wing
    
    % passing the each row of the wake to the next level
    for i=1:tmax-1
        Gw(1+2*i*Ny:2*(i+1)*Ny)=Gw(1+2*(i-1)*Ny:2*i*Ny) ;
    end
    
    Gw(1:2*Ny)= G(1+2*(Nx-1)*Ny:2*Nx*Ny,t);
    
    [cLnewv(t),cl_section(:,t),dp(:,:,t)]=force_calc(Nx,Ny,a,G,Gs(:,1),t,dt,S,dl_x,dly,uinf,aoa,Lam,dih,aoaf_L,aoaf_R);  
end
toc
%% figure
%plot(linspace(0,vtmax*dt*uinf/cr,vtmax),cLnewv);xlabel('time'); ylabel('C_L');grid minor;hold on;
%figure
%plot(y,cl_section(:,10));
%axis([0 9 0 0.55]);
% surf(x,y,z,[G(1:2*Ny,34)';G(2*Ny+1:4*Ny,34)';G(4*Ny+1:6*Ny,34)';G(6*Ny+1:8*Ny,34)'])
figure
surf(repmat(linspace(0,vtmax*dt*uinf/cr,vtmax),2*Ny,1),repmat(y(1,1:2*Ny),vtmax,1)',cl_section)
xlabel('time');ylabel('span');zlabel('c_l');
view(-150,15);

%% plotting real wing -- converting the middle points data to edges.
bn=b+0.5*b/Ny/2;
ARn=bn^2/(S*(Nx+1)/Nx*(Ny+1/4)/Ny);
trn=tr;
[xi,yi,zi,xcoli,ycoli,zcoli,ni,dl_xi,dlyi,Si,alphai,cri]=geometry(ARn,bn,trn,Nx+1,2*Ny+1/2,Lam,dih,aoa,0,0);

xi(Nx+1,4*Ny-2)=NaN;xi(Nx+1,4)=NaN;
%xi(Nx,4*Ny-2)=NaN;

% right flap geometry
xifr(1,1)=xi(Nx,3);xifr(1,2)=xi(Nx,5);
xifr(2,1)=xi(Nx,3)+dl_xi(3)*cos(aoa+aoaf_R);
xifr(2,2)=xi(Nx,5)+dl_xi(5)*cos(aoa+aoaf_R);
yifr(1,1)=yi(Nx,3);yifr(1,2)=yi(Nx,5);yifr(2,1)=yi(Nx,3);yifr(2,2)=yi(Nx,5);
zifr(1,1)=zi(Nx,3);zifr(1,2)=zi(Nx,5);
zifr(2,1)=zi(Nx,3)-dl_xi(3)*sin(aoaf_R+aoa)*cos(dih);
zifr(2,2)=zi(Nx,5)-dl_xi(5)*sin(aoaf_R+aoa)*cos(dih);

% left flap geometry
xifl(1,1)=xi(Nx,4*Ny-3);xifl(1,2)=xi(Nx,4*Ny-1);
xifl(2,1)=xi(Nx,4*Ny-3)+dl_xi(4*Ny-3)*cos(aoa+aoaf_L);
xifl(2,2)=xi(Nx,4*Ny-1)+dl_xi(4*Ny-1)*cos(aoa+aoaf_L);
yifl(1,1)=yi(Nx,4*Ny-3);yifl(1,2)=yi(Nx,4*Ny-1);yifl(2,1)=yi(Nx,4*Ny-3);yifl(2,2)=yi(Nx,4*Ny-1);
zifl(1,1)=zi(Nx,4*Ny-3);zifl(1,2)=zi(Nx,4*Ny-1);
zifl(2,1)=zi(Nx,4*Ny-3)-dl_xi(4*Ny-3)*sin(aoaf_L+aoa)*cos(dih);
zifl(2,2)=zi(Nx,4*Ny-1)-dl_xi(4*Ny-1)*sin(aoaf_L+aoa)*cos(dih);

%%
for i=1:Nx
    
    xv(1+(i-1)*2*Ny:2*Ny*i)=x(i,1:2*Ny);
    yv(1+(i-1)*2*Ny:2*Ny*i)=y(i,1:2*Ny);
    zv(1+(i-1)*2*Ny:2*Ny*i)=z(i,1:2*Ny);
    for j=1:vtmax
        dpv(1+(i-1)*2*Ny:2*Ny*i,j)=dp(i,1:2*Ny,j);
    end
end
% the first is the interpolatoin (linear) and the second is extraploation (nearest)
aa=scatteredInterpolant(xv',yv',zv',dpv(:,5),'linear','nearest')
figure;
surf(xi(:,1:2*Ny+1),yi(:,1:2*Ny+1),zi(:,1:2*Ny+1),aa(xi(:,1:2*Ny+1),yi(:,1:2*Ny+1),zi(:,1:2*Ny+1)));hold on;axis('equal');
surf( fliplr(xi(:,2*Ny+1:4*Ny+1)), fliplr(yi(:,2*Ny+1:4*Ny+1)), fliplr(zi(:,2*Ny+1:4*Ny+1)),aa( fliplr(xi(:,2*Ny+1:4*Ny+1)), fliplr(yi(:,2*Ny+1:4*Ny+1)), fliplr(zi(:,2*Ny+1:4*Ny+1))))

surf(xifr,yifr,zifr,aa(xifr,yifr,zifr));surf(xifl,yifl,zifl,aa(xifr,yifr,zifr));

%surf(xifr,yifr,zifr,repmat(dp(Nx,2,10),2,2));surf(xifl,yifl,zifl,repmat(dp(Nx,2*Ny-1,10),2,2));

set(gca, 'CLim', [0, .9]);
colorbar;
% figure;
% surf(x,y,z,dp(:,:,20));

%% movie

% figure
% for i=1:vtmax
%   aa=scatteredInterpolant(xv',yv',zv',dpv(:,i));
% surf(xi,yi,zi,aa(xi,yi,zi));hold on;axis('equal');
% surf(xifr,yifr,zifr,aa(xifr,yifr,zifr));surf(xifl,yifl,zifl,aa(xifr,yifr,zifr));
%
% %     surf(x,y,z,dp(:,:,i));drawnow;
%    set(gca, 'CLim', [0, .5]);colorbar;
%      axis('equal');
%     Mv(i)=getframe;
% end
% v = VideoWriter('movie.avi');
% open(v)
% writeVideo(v,Mv)
% close(v)
