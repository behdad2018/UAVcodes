% simplified unsteady 3D wing lifting line method by B. Davoudi
% Aerospace Engineering Department, University of Michigan 7/27/2016

tic
clc; clear all; close all;
Nx=4;                        % chordwise panel per half wing
Ny=12;                        % spanwise panel per hald wing
AR=8;                        % Aspect ratio
b=1;                         % wing span
tr=1;                      % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoa=degtorad(8);             % angle of attack
uinf=25;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector
aoaf_L=degtorad(10);          % flap angle left positive y looking fron front to airplane
aoaf_R=degtorad(0);         % flap angle right negative y looking fron front to airplane

% geometry
[x,y,z,xcol,ycol,zcol,n,dl_x,dly,S,alpha,cr]=geometry(AR,b,tr,Nx,Ny,Lam,dih,aoa,aoaf_L,aoaf_R);
                                                    
 figure(1);
 surf(x,y,z);axis('equal')

% solving for vorticty
[Gs,Am_wing,~,a,a_d,w_ind_drag]=fast_steady(x,y,z,xcol,zcol,n,dl_x,dly,Nx,0,Ny,0,u,alpha,Lam,dih,b);

% calculate forces
[cLnewvs,cl_section,cd_section,dps,d,roll,yaw]=force_calc(Nx,Ny,a,a_d,w_ind_drag,b,Gs,Gs,1,1,S,dl_x,dly,uinf,aoa,Lam,dih,aoaf_L,aoaf_R);
                                 %   (Nx,Ny,a,G,G5,t,dt,S,dl_x,dly,uinf,aoa,Lam,dih,aoaf_L,aoaf_R)
figure(2)
plot(y,cl_section) 
xlabel('y (wing span)');ylabel('C_l');hold on

figure(3)
plot(y,cd_section) 
xlabel('y (wing span)');ylabel('C_d');hold on

figure(6)
plot(y,cl_section./cd_section) 
xlabel('y (wing span)');ylabel('Cl/cd');hold on

%%  plotting real wing -- converting the middle points data to edges.
%   zr=bw/16+yr1(:,Nyr+1:end);xr=bw/2+xr1(:,Nyr+1:end);yr=zr1(:,Nyr+1:end);

figure(4)

sp(dps,x,y,z,b,b,Nx,Ny,tr,S,Lam,dih,aoa,aoaf_L,aoaf_R,1,max(max(dps)),0,0);

toc

figure(5)

sp(d,x,y,z,b,b,Nx,Ny,tr,S,Lam,dih,aoa,aoaf_L,aoaf_R,1,max(max(d)),0,0);