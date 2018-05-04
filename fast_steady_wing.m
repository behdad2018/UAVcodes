function [G,Am_wing,a]=fast_steady_wing(x,y,z,xcol,zcol,n,dl_x,dly,Nx,Ny,u,alpha,Lam,dih,b)
% genrating the wake for a steady flight
% finding the steady solution with the panel method
% original code to do only wing


Awing=zeros(2*Nx*Ny,2*Nx*Ny);
Awake=zeros(2*Nx*Ny,2*Nx*Ny);

for i=1:2*Nx*Ny    
    index1col=ceil(i/(2*Ny));
    index2col=i-(2*Ny)*(index1col-1);
    xcol1=xcol(index1col,index2col);
    ycol1=y(index1col,index2col);
    zcol1=zcol(index1col,index2col);
    %  for wing
    for j=1:2*(Nx+1)*Ny  % panels interations on the wing -- finding the matrix A
        index1=ceil(j/(2*Ny));
        index2=j-(2*Ny)*(index1-1);
        % for wing
        if  index1<=Nx
            x1=x(index1,index2);
            y1=y(index1,index2);
            z1=z(index1,index2);
            [~,Awing(i,j)]=vortexring(n(:,i),dl_x(index2),dly,alpha(index1,index2),Lam,dih,xcol1,ycol1,zcol1,x1,y1,z1,1);       % Effect of ring vortices on the wing
        % for wake
        elseif  index1>Nx
            xw1=x(end,index2)+dl_x(index2)*cos(alpha(index1-1,index2));
            yw1=y(end,index2);
            zw1=z(end,index2)-dl_x(index2)*sin(alpha(index1-1,index2))*cos(dih);
            [a1(i,j-2*Ny,:),Awake(i,j-2*Ny)]=vortexring(n(:,i),20*b,dly,0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,1);
        end
    end
end
A=Awing+Awake;
% the inverse matrix A that will be used in the unsteady solution
Am_wing=inv(Awing);
RHS=dot(-repmat(u',1,2*Nx*Ny),n)';
G=A\RHS;    % vorticity on the wing
%these are the induced velocty due to the wake on the wing panels
a(1,:)=a1(:,:,1)*G;
a(2,:)=a1(:,:,3)*G;
a(3,:)=a1(:,:,3)*G;
end
