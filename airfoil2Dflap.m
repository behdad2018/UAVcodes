% lift prediction using 2-D airfoilf theory
c=1;                % chord
xf=5/6;            % flap location, x is non-dimensional
aoa=degtorad(5);    % angle of attack
% aoaf             flap angle of attack

% simplified unsteady 3D wing lifting line method by B. Davoudi 7/27/2016

Nx=6;                        % chordwise panel per hald wing
Ny=50;                        % spanwise panel per hald wing
AR=100000;                        % Aspect ratio
b=1;                         % wing span
tr=1;                        % taper ratio
Lam=degtorad(0);             % sweep angle, backward swept, positive
dih=degtorad(0);             % dihedral angle defined at the c/4
aoa=degtorad(5);             % angle of attack
uinf=1;                      % incidence velocity
u=uinf*[1 0 0];              % incidence velocity vector


for i=1:31

aoaf=degtorad(i-1)   % angle of attack

% trasnformation
thf=acos(1-2*xf);

% dz/dx of the flap...note dz/dx of the rest is zero.

th=linspace(thf,pi,100);
dzdx=repmat(-tan(aoaf),1,100);

a0=aoa-(1/pi)*trapz(th,dzdx);
a1=2/pi*trapz(th,dzdx.*cos(th));


cl2D(i)=2*pi*(a0+0.5*a1);


aoaf_L=aoaf;
aoaf_R=aoaf;

% geometry
[x,y,z,xcol,ycol,zcol,n,dl_x,dly,S,alpha,cr]=geometry(AR,b,tr,Nx,Ny,Lam,dih,aoa,aoaf_L,aoaf_R);


[G15,Awake,a]=fast_steady(x,y,z,xcol,zcol,n,dl_x,dly,Nx,Ny,u,alpha,Lam,dih,b);

[cLnewvs,cls_section(:,i),dps]=force_calc(Nx,Ny,a,G15,G15,1,1,S,dl_x,dly,uinf,aoa,Lam,dih,aoaf_L,aoaf_R);


end
plot([0:30],cl2D);hold on;plot([0:30],cls_section(50,:));grid on;
xlabel('\alpha_f');ylabel('c_l');
legend('2D airfoil theory','3D vortex method, AR=1000');
