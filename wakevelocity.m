% This function obtained wake velocity of all of the wake pannels
% by considering the induced velocity by wake pannels and wing pannels

function [u]=wakevelocity(Nx,Ny,Nw,dlw,dlx,dly,aoa,Lam,dih,x,y,z,xw,yw,zw,G,Gw)

for i=1:2*Nw*Ny
    
    % coordinates of the wake pannel the convection speed is to be found
    ind1=ceil(i/(2*Ny));
    ind2=i-(2*Ny)*(ind1-1);
    
    u1=0;
    
    % the induced velcoity is found at the center of the vortex panel]
    % note that the apnnels remain flat wrt free stream and won't rotate
    dlw1=dlw(ind1,ind2);
    xw1=xw(ind1,ind2) + 0.5*dlw1;
    yw1=yw(ind1,ind2);
    zw1=zw(ind1,ind2);
    
    % induced velocity by other wake pannels
    for k=1:2*Nw*Ny
        
            in1=ceil(k/(2*Ny));
            in2=k-(2*Ny)*(in1-1);
            
            xw2=xw(in1,in2);
            yw2=yw(in1,in2);
            zw2=zw(in1,in2);
    
            [u2,a]=vortexring([0 0 1],dlw(in1,in2),dly,0,Lam,dih,xw1,yw1,zw1,xw2,yw2,zw2,Gw(k));
            u1=u1+u2;
        
    end
    
    % induced velocity by wing pannels
    for k=1:2*Nx*Ny
        
        in1=ceil(k/(2*Ny));
        in2=k-(2*Ny)*(in1-1);
        
        x1=x(in1,in2);  
        y1=y(in1,in2); 
        z1=z(in1,in2);
        
        [u2,a]=vortexring([0 0 1],dlx,dly,aoa,Lam,dih,xw1,yw1,zw1,x1,y1,z1,G(k));
        
        u1=u1+u2;
        
    end
    
    u(:,ind1,ind2)=u1;
    
end

end