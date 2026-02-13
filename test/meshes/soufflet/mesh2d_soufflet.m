% Channel: regular mesh, equilateral triangles
% Soufflet et al. 20,10,5 and 2 km meshes.

% Cyclic_length=4.5
% ============
% Define the number of levels
nl=41;                 % 41 for 10 km, 61 for 5 km
H=4000;                % Depth in m
alpha=1.1;            % The next layer is alpha times thicker
                       % 1.1 for 10 km, 1.06 for 5 km    
% The layer thicknesses will be computed based on this information.                

% Domain size

re=6400;               % Earth's radius
Lx=(4.5-2*0.09)*pi*re/180;      % In km, add 10km in deg so I add one longitude column
                       % We ensure that the cyclic_length is 4.5, othervise  
                       % it is very close to 500 km of Soufflet et al.
Ly=2010-2*8.75;	       % adding (lat2-lat1)*pi*re/180 = 8.75 (height of one triangle in km) for a 10 km mesh
  dx=20;               % Approximate mesh resolution. Change it here.
  nn=floor(Lx/dx);
  dx=Lx/nn;            % dx is close to 10 km, and we have the integer number
                       % of intervals within Lx
  
  dy=dx*sqrt(3)/2;
  nn=floor(Ly/dy);
  dy=Ly/nn;            % we ensure that the channel width is exactly Ly
   
  dx=dx*180/re/pi;     % In degrees    
  dy=dy*180/re/pi;
  
  
  lon=0:dx:(4.5-2*0.09);
  lat=0:dy:180*Ly/re/pi;

  disp(lon);
  disp(lat);
  
  nx=length(lon);
  ny=length(lat);
  nodnum=1:nx*ny;
  nodnum=reshape(nodnum,[ny, nx]);
  xcoord=zeros([ny, nx]);
  xnum=zeros([ny, nx]);
  ycoord=xcoord;
  ynum=zeros([ny, nx]);
  for n=1:nx,
  ycoord(:,n)=lat';
  ynum(:,n)=(1:ny)';
  end;
   
  for n=1:ny,
  xcoord(n,:)=lon;
  xnum(n,:)=(1:nx);
  end;
  for n=2:2:ny,
  xcoord(n,:)=xcoord(n,:)+0.5*dx;
  end;
  xcoord=reshape(xcoord,[1,nx*ny]);
  ycoord=reshape(ycoord,[1,nx*ny]); 
  
  tri=[];
  for n=1:nx-1,
      for nn=1:2:ny-1
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn,n+1)]];
      tri=[tri; [nodnum(nn+1,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
      end;
      for nn=2:2:ny-1
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n),nodnum(nn+1,n+1)]];
      tri=[tri; [nodnum(nn,n),nodnum(nn+1,n+1),nodnum(nn,n+1)]];
      end;
  end;    
  % plot
  trimesh(tri, xcoord', ycoord', zeros(nx*ny,1)); 
  view(2)
   
  % Cyclic reduction:
  % The last ny nodes in xcoord, ycoord are equivalent to 
  % the first ny nodes.
  ai=find(tri>(nx-1)*ny);
  tri(ai)=tri(ai)-(nx-1)*ny;
  xcoord=xcoord(1:(nx-1)*ny);
  ycoord=ycoord(1:(nx-1)*ny);
 
  % Cyclic reduction means that the last column has to be removed 
  % in xnum, ynum 
  xnum=xnum(:, 1:nx-1);
  ynum=ynum(:, 1:nx-1);
  xnum=reshape(xnum,[ny*(nx-1),1]);
  ynum=reshape(ynum,[ny*(nx-1),1]);
  
  n2d=(nx-1)*ny;
  nodes=zeros([4, n2d]);
  nodes(1,:)=1:n2d;
  nodes(2,:)=xcoord;
  nodes(3,:)=ycoord;
  nodes(4,:)=zeros(size(ycoord));
  % Set indices to 1 on vertical walls
  ai=find(ycoord==min(lat));
  nodes(4,ai)=1;
  ai=find(ycoord==max(lat));
  nodes(4,ai)=1;
  
  % Define levels:

  dz=H*(1-alpha)/(1-alpha^(nl-1))       % dz of the top layer; then alpha*dz,...
zbar=zeros([1,nl]);
zbar(2)=dz; 
for n=3:nl
zbar(n)=zbar(n-1)*alpha;
end
for n=1:nl-1,
    zbar(n+1)=zbar(n)+zbar(n+1);
end
dd=-H*ones(size(xcoord));

  % Output 2D mesh 
  fid = fopen('/albedo/work/projects/p_clidyn_work/rjuhrban/channel/m20/nod2d.out','w');
        fprintf(fid,'%8i \n',n2d);
        fprintf(fid,'%8i %8.4f %8.4f %8i\n',nodes); 
        fclose(fid);
 clear nodes       
        
 fid=fopen('/albedo/work/projects/p_clidyn_work/rjuhrban/channel/m20/elem2d.out','w');
        fprintf(fid,'%8i \n', length(tri(:,1)));
        fprintf(fid,'%8i %8i %8i\n',tri');
        fclose(fid);
 
 
      fid=fopen('/albedo/work/projects/p_clidyn_work/rjuhrban/channel/m20/aux3d.out', 'w');
      fprintf(fid,'%g\n', nl);
      fprintf(fid,'%g\n', zbar);
      fprintf(fid,'%7.1f\n', dd);
       
fclose(fid);
