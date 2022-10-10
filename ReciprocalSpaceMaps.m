function ReciprocalSpaceMaps(wl, twoth_d, phi_d, D, p, Ny, Nx, n0, m0, nb, mb, folder, fnamepng, sclfactor, Az, El, th0, dth)

tifFiles = dir([folder '\*.tif']);
N = length(tifFiles);    % number of image files
th =  th0+(0:N-1)*dth;   % array of th values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Common matrices to all th values          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sd = [cosd(twoth_d)*[cosd(phi_d) sind(phi_d)] sind(twoth_d)];
xd = [sind(phi_d) -cosd(phi_d) 0];
yd = [-sind(twoth_d)*[cosd(phi_d) sind(phi_d)] cosd(twoth_d)];
ry = -([1:Ny] - n0)'*p;  % vertical pixels, 1st at the upper left corner
rx = ([1:Nx] - m0)*p;    % horizontal pixels
Rx = D*sd(1) + repmat(ry*yd(1),1,Nx) + repmat(rx*xd(1),Ny,1);
Ry = D*sd(2) + repmat(ry*yd(2),1,Nx) + repmat(rx*xd(2),Ny,1);
Rz = D*sd(3) + repmat(ry*yd(3),1,Nx) + repmat(rx*xd(3),Ny,1);
R = sqrt(Rx.*Rx + Ry.*Ry + Rz.*Rz);
aux = (2*pi/wl);
 Qx = aux*(Rx./R - 1);
 Qy = aux*Ry./R;
 Qz = aux*Rz./R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D matrices                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DQx = zeros(Ny,Nx,N); DQy = DQx; DQz = DQx; A = DQx;
for j=1:N
  DQx(:,:,j) = Qx.*repmat(cosd(th(j)),Ny,Nx) + Qz.*repmat(sind(th(j)),Ny,Nx);  
  DQy(:,:,j) = Qy;
  DQz(:,:,j) = Qx.*repmat(-sind(th(j)),Ny,Nx) + Qz.*repmat(cosd(th(j)),Ny,Nx);
  Ij = double(imread([folder '\' tifFiles(j).name])); Ij(nb,mb) = 1;
  A(:,:,j) = Ij;
end

A(A<1)=1;
Amax=max(max(max(A)));
whalf = log10(0.5*Amax);
Z = log10(A);
hf1 = figure(1);
clf
set(hf1,'InvertHardcopy','off','Color','w')
patch(isosurface(DQx,DQy,DQz,Z,sclfactor(1)*whalf),'FaceColor','red',...
           'EdgeColor','none','FaceAlpha',1);
hold on
patch(isosurface(DQx,DQy,DQz,Z,sclfactor(2)*whalf),'FaceColor','green',...
           'EdgeColor','none','FaceAlpha',.22);
patch(isosurface(DQx,DQy,DQz,Z,sclfactor(3)*whalf),'FaceColor','blue',...
           'EdgeColor','none','FaceAlpha',.04);
hold off
       
view(3)
daspect([1 1 1])
axis tight
camup([0 0 1])
camlight
lighting phong
set(gca,'Color',[0.97 0.97 0.97],'Box','on',...
    'FontSize',10,'LineWidth',1,'FontName','Arial')
xlabel(['\Delta Q_x (' char(197) '^{-1})'],'FontSize',12)
ylabel(['\Delta Q_y (' char(197) '^{-1})'],'FontSize',12)
zlabel(['\Delta Q_z (' char(197) '^{-1})'],'FontSize',12)

camlight
axis tight
view(Az,El); 
grid; 
set(gca,'XLim',[-.05 .05],'YLim',[-.05 .05])
print('-dpng','-r300',[fnamepng '.png'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D projection                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qz x Qx 
hf2 = figure(2); clf
set(hf2,'InvertHardcopy','off','Color','w')
XX=zeros(Ny,N);
ZZ=XX;
II=XX;
for n=1:N
    S = zeros(Ny,1);
    for m=1:Nx
        S = S + A(:,m,n);
    end
      II(:,n) = S; 
end
IImax = 0.5*max(max(II));
II(II>IImax) = IImax;
for n=1:N, XX(:,n) = DQx(:,m0,n); end
for n=1:N, ZZ(:,n) = DQz(:,m0,n); end
contourf(XX,ZZ,log10(II),12)
set(gca,'FontName','Arial','FontSize',10,'LineWidth',1);
axis image
axis tight
c = colorbar;
c.Label.String = 'log(counts)';
xlabel(['\Delta Q_x (' char(197) '^{-1})'],'FontSize',12)
ylabel(['\Delta Q_z (' char(197) '^{-1})'],'FontSize',12)
print('-dpng','-r300',[fnamepng 'QxQz.png'])
