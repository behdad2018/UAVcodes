% simplified unsteady 3D wing lifting line method by B. Davoudi
% Aerospace Engineering Department, University of Michigan 7/27/2016
% updated to incorporate the rudder 9/7/2017
function [G,Am_wing,A,a,a_d,w_ind_drag]=fast_steady_tail_rudder(x,y,z,xcol,ycol,zcol,n,dl_x,dly,Nxw,Nxt,Nxr,Nyw,Nyt,Nyr,u,alpha,Lam,dih,b)
% genrating the wake for a steady flight
% finding the steady solution with the panel method

Awing=zeros(2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr,2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr);
Awake=zeros(2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr,2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr);
a1=zeros(2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr,2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr,3);

%for i=1:2*Nx*Ny

for i=1:2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr  % total of the panels in the wing and wake
    
    % once the panels on the wing are considered, tail panels are the next
    % starting from 2*Nxw*Nyw+1 to 2*Nxw*Nyw+2*Nxt*Nyt
    
    % for wing
    if i<=2*Nxw*Nyw
        index1col=ceil(i/(2*Nyw));
        index2col=i-(2*Nyw)*(index1col-1);      
%        once we start consdiering the tail panels
    elseif i>2*Nxw*Nyw && i<=2*Nxw*Nyw+2*Nxt*Nyt
        index1col=Nxw+ceil((i-2*Nxw*Nyw)/(2*Nyt));                     % goes from Nxw+1 to Nxt
        index2col=i-2*Nxw*Nyw-(2*Nyt)*(index1col-Nxw-1);               % goes from 1 to Nyt
        % once we start consdiering the rudder panels
    elseif i>2*Nxw*Nyw+2*Nxt*Nyt
        index1col=Nxw+Nxt+ceil((i-2*Nxw*Nyw-2*Nxt*Nyt)/Nyr);       % goes from Nxr+1 to Nxr
        index2col=i-2*Nxw*Nyw-2*Nxt*Nyt-Nyr*(index1col-Nxw-Nxt-1);      % goes from 1 to Nyr
    end

    % collocation points for both wing and tail
    xcol1=xcol(index1col,index2col);
    ycol1=ycol(index1col,index2col);
    zcol1=zcol(index1col,index2col);
    
    % correcting the position of the collocation points on rudder
    if i>2*Nxw*Nyw+2*Nxt*Nyt 
     M=[1,0,0;0,cosd(90),-sind(90);0,sind(90),cosd(90)]*[xcol1;ycol1;zcol1];
     xcol1=M(1);ycol1=M(2);zcol1=M(3);     
    end
    
    %% for finding matrix A only panels on the wing, tail and rudder
    for j=1:2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr
        
    if j<=2*Nxw*Nyw
        index1=ceil(j/(2*Nyw));
        index2=j-(2*Nyw)*(index1-1); 
        index3=index2; 
        theta=0;
        onetwothree=1;
        % once we start consdiering the tail panels
    elseif j>2*Nxw*Nyw && j<=2*Nxw*Nyw+2*Nxt*Nyt
        index1=Nxw+ceil((j-2*Nxw*Nyw)/(2*Nyt));                     % goes from Nxw+1 to Nxt
        index2=j-2*Nxw*Nyw-(2*Nyt)*(index1-Nxw-1);               % goes from 1 to Nyt
        index3=index2+2*Nyw;
        theta=0;
        onetwothree=2;
        % once we start consdiering the rudder panels
    elseif j>2*Nxw*Nyw+2*Nxt*Nyt
        index1=Nxw+Nxt+ceil((j-2*Nxw*Nyw-2*Nxt*Nyt)/Nyr);       % goes from Nxt+1 to Nxr
        index2=j-2*Nxw*Nyw-2*Nxt*Nyt-Nyr*(index1-Nxw-Nxt-1);      % goes from 1 to Nyr
        index3=index2+2*Nyw+2*Nyt;
        theta=90;
        dih=0;
        onetwothree=3;
    end
 
        x1=x(index1,index2);
        y1=y(index1,index2);
        z1=z(index1,index2);

        [~,Awing(i,j)]=vortexring(n(:,i),dl_x(index3),dly(onetwothree),alpha(index1,index2),Lam,dih,xcol1,ycol1,zcol1,x1,y1,z1,1,theta);       % Effect of ring vortices on the wing
        % for drag calculation
        [~,Awing_drag(i,j)]=vortexring_drag(n(:,i),dl_x(index3),dly(onetwothree),alpha(index1,index2),Lam,dih,xcol1,ycol1,zcol1,x1,y1,z1,1,theta);
    end
    %% wake effect to be included in the A matrix
    % Note: we have two wakes one from the wing and one from the tail. Also
    % we think of the wake as if it were extended of the last panel
    % wing:
    for j=1:2*Nyw
        
        % first row of the wake, index1=Nxw
        xw1=x(Nxw,j)+dl_x(j)*cos(alpha(Nxw,j));
        yw1=y(Nxw,j);
        zw1=z(Nxw,j)-dl_x(j)*sin(alpha(Nxw,j))*cos(dih);
        theta=0;
        [a1(i,2*(Nxw-1)*Nyw+j,:),Awake(i,2*(Nxw-1)*Nyw+j)]=vortexring(n(:,i),20*b,dly(1),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,theta);
        [a1_d(i,2*(Nxw-1)*Nyw+j,:),~]=vortexring_drag(n(:,i),20*b,dly(1),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,theta);
    end
    % tail:
    for j=1:2*Nyt
        
        % first row of the wake index1=Nxt
        xw1=x(Nxt+Nxw,j)+dl_x(2*Nyw+j)*cos(alpha(Nxt+Nxw,j));
        yw1=y(Nxt+Nxw,j);
        zw1=z(Nxt+Nxw,j)-dl_x(2*Nyw+j)*sin(alpha(Nxt+Nxw,j))*cos(dih);
        theta=0;
        [a1(i,2*Nxw*Nyw+2*(Nxt-1)*Nyt+j,:),Awake(i,2*Nxw*Nyw+2*(Nxt-1)*Nyt+j)]=vortexring(n(:,i),20*b,dly(2),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,theta);
        [a1_d(i,2*Nxw*Nyw+2*(Nxt-1)*Nyt+j,:),~]=vortexring(n(:,i),20*b,dly(2),0,Lam,dih,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,theta);
    end
    
        % rudder:
    for j=1:Nyr

        xw1=x(Nxt+Nxw+Nxr,j)+dl_x(2*Nyw+2*Nyt+j)*cos(alpha(Nxt+Nxw+Nxr,j));
        yw1=y(Nxt+Nxw+Nxr,j);
        zw1=z(Nxt+Nxw+Nxr,j)-dl_x(2*Nyw+2*Nyt+j)*sin(alpha(Nxt+Nxw+Nxr,j))*cos(0);
        theta=90;
       [a1(i,2*Nxw*Nyw+2*Nxt*Nyt+(Nxr-1)*Nyr+j,:),Awake(i,2*Nxw*Nyw+2*Nxt*Nyt+(Nxr-1)*Nyr+j)]=vortexring(n(:,i),20*b,dly(3),0,0,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,theta);
       [a1_d(i,2*Nxw*Nyw+2*Nxt*Nyt+(Nxr-1)*Nyr+j,:),~]=vortexring(n(:,i),20*b,dly(3),0,Lam,0,xcol1,ycol1,zcol1,xw1,yw1,zw1,1,theta);
    end
    
end

A=Awing+Awake;

% the inverse matrix A that will be used in the unsteady solution

Am_wing=inv(Awing);

U=repmat(u',1,2*Nxw*Nyw+2*Nxt*Nyt+Nxr*Nyr);

%%

RHS=dot(-U,n)';

G=A\RHS;    % vorticity on the wing
size(G);

% these are the induced velocty due to the wake on the wing panels
a(1,:)=a1(:,:,1)*G;
a(2,:)=a1(:,:,2)*G;
a(3,:)=a1(:,:,3)*G;


a_d(1,:)=a1_d(:,:,1)*G;
a_d(2,:)=a1_d(:,:,2)*G;
a_d(3,:)=a1_d(:,:,3)*G;

w_ind_drag=Awing_drag*G;

end



