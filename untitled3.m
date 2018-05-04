% location of lifting line from -y to y
%leng=linspace(0.5*dly,0.5*wl-0.5*dly,Ny);

% for j=1:2*Ny
%     leng(j)=-cos((j-1)*pi/(2*Ny-1))*0.5*b;
%     if j<2*Ny
%         lengcol(j)=-cos((j-0.5)*pi/(2*Ny-1))*0.5*b;
%     end
% end
% 
% dl_x=c/Nx;                                        % pannel length stream wise
% dlx=min(dl_x);                                    % this should be modified for tr~=1
% 
% % locations of lifting line on the wing (Nx=1)
% for i=1:Nx
%     y(i,:)=leng*cos(Lam)*cos(dih);
%     z(i,:)=leng*(sin(dih)-sin(aoa)*sin(Lam))-dlx*(i-1)*sin(aoa);
%     x(i,:)=leng*sin(Lam)*cos(aoa) + dlx*(i-1)*cos(aoa);
%     
%     ycol(i,:)=lengcol*cos(Lam)*cos(dih);
%     zcol(i,:)=lengcol*(sin(dih)-sin(aoa)*sin(Lam))-dlx*(i-1)*sin(aoa);
%     xcol(i,:)=lengcol*sin(Lam)*cos(aoa) + dlx*(i-1)*cos(aoa);
% end
% 
% dly=diff(y);               % pannel length span wise