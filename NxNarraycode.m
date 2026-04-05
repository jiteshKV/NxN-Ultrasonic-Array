% Simulations In Three Dimensions for NxN array

clearvars; 

% ========= Grid Set Up ===================================================

% create the computational grid
Nx = 65;            % number of grid points in the x direction
Ny = 65;            % number of grid points in the y direction
Nz = 65;            % number of grid points in the z direction
dx = 0.1e-2;        % grid point spacing in the x direction [m]
dy = 0.1e-2;        % grid point spacing in the y direction [m]
dz = 0.1e-2;        % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% Focus point coordinates in meters
fX=0.024;%0.032;
fY=-0.024;%0.015;
fZ=-0.024;%0.015;
narray=5;   % N in NxN arrat




% ========= Medium Set Up ===================================================
% define the properties of the propagation medium
medium.sound_speed = 343; %* ones(Nx, Ny, Nz);	% [m/s]


% create the time array
kgrid.makeTime(medium.sound_speed);





% ========= Source Set Up ===================================================

% -----define a NxN Array ------
source1.p_mask = zeros(Nx, Ny,Nz); % Mask for 1st NxN array
source2.p_mask=source1.p_mask; % Mask for 2nd NxN array
source3.p_mask=source1.p_mask; % Mask for both NxN array combined
gapEl=3;


%----------Combined Array -----------------
midElZ=ceil(Nz/4); 
midElY=ceil(Ny/2);
spreadEl=(1:narray)*gapEl;
spreadEl=spreadEl-spreadEl(1,ceil(size(spreadEl,2)/2));
rowPos=spreadEl+midElZ;
colPos=spreadEl+midElY;
source3.p_mask(2,colPos,rowPos)=1;

midElZ=ceil(3*Nz/4);
midElY=ceil(Ny/2);
spreadEl=(1:narray)*gapEl;
spreadEl=spreadEl-spreadEl(1,ceil(size(spreadEl,2)/2));
rowPos=spreadEl+midElZ;
colPos=spreadEl+midElY;
source3.p_mask(2,colPos,rowPos)=1;

%Determining Time delay
[cart_data, order_index] = grid2cart(kgrid, source3.p_mask);
focalPt=repmat([fX;fY;fZ],1,size(cart_data,2));
dist=sqrt(sum((focalPt-cart_data).^2));
pathDiff=dist-min(dist);
timeDelay=pathDiff./medium.sound_speed;
timeDelay=timeDelay';
timeDelay=repmat(timeDelay,1,kgrid.Nt);

% define a time varying sinusoidal source
source_freq = 4e4;  % [Hz]
source_mag = 1;     % [Pa]
tArray=repmat(kgrid.t_array,size(cart_data,2),1); 
source3.p = source_mag * sin(2 * pi * source_freq * (tArray+timeDelay));

%-------Setting Bottom Array--------
midElZ=ceil(Nz/4); 
midElY=ceil(Ny/2);
spreadEl=(1:narray)*gapEl;
spreadEl=spreadEl-spreadEl(1,ceil(size(spreadEl,2)/2));
rowPos=spreadEl+midElZ;
colPos=spreadEl+midElY;
source1.p_mask(2,colPos,rowPos)=1;

% define a time varying sinusoidal source
source_freq = 4e4;  % [Hz]
source_mag = 1;     % [Pa]
tArray=repmat(kgrid.t_array,((narray)^2),1); 
td_bot=timeDelay((1:(narray)^2),:);
source1.p = source_mag * sin(2 * pi * source_freq * (tArray+td_bot));

%----------Setting Top Array-----------------

midElZ=ceil(3*Nz/4);
midElY=ceil(Ny/2);
spreadEl=(1:narray)*gapEl;
spreadEl=spreadEl-spreadEl(1,ceil(size(spreadEl,2)/2));
rowPos=spreadEl+midElZ;
colPos=spreadEl+midElY;
source2.p_mask(2,colPos,rowPos)=1;

% define a time varying sinusoidal source
source_freq = 4e4;  % [Hz]
source_mag = 1;     % [Pa]
tArray=repmat(kgrid.t_array,((narray)^2),1); 
source2.p = source_mag * sin(2 * pi * source_freq * (tArray+timeDelay((narray)^2+1:end,:)));




% ========= Sensor Set Up ===================================================

% Defining Sensor 1 - Axial

sensor1.mask = zeros(Nx, Ny,Nz);
sensor1.mask(:,midElY,:)=1;
sensor1.record = {'p', 'p_rms'};

% Defining Sensor 2 - Transverse
[fpoint_mask, order_index, reorder_index] = cart2grid(kgrid, [fX;fY;fZ]);
[M,fidx]=max(fpoint_mask(:));
[row,col,page]=ind2sub([Nx Ny Nz],fidx);
sensor2.mask = zeros(Nx, Ny,Nz);
sensor2.mask(row,:,:)=1;
sensor2.record = {'p', 'p_rms'};



% ========= Simulation Set Up ===================================================

% input arguments
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'PlotSim', false, ...
    'DisplayMask', 'off'};
% run the simulation
sensor_data11 = kspaceFirstOrder3D(kgrid, medium, source1, sensor1, input_args{:});
sensor_data21 = kspaceFirstOrder3D(kgrid, medium, source2, sensor1, input_args{:});
sensor_data31 = kspaceFirstOrder3D(kgrid, medium, source3, sensor1, input_args{:});
sensor_data32 = kspaceFirstOrder3D(kgrid, medium, source3, sensor2, input_args{:});




% =========================================================================
% VISUALISATION
% =========================================================================


fposxz=[fX-dx,fZ-dz,2*dx,2*dz]; % focal position rectangle at XZ sensor plane
fposyz=[fY-dy,fZ-dz,2*dy,2*dz]; % focal position rectangle at YZ sensor plane


%--------------Top Array Only --------------------
shg;figure('Position', [187,349,1422,558]);
sgt=sgtitle('Top Array Only');sgt.FontSize=20;
%plot 1
subplot(1,3,1);voxelPlot(double(source2.p_mask | sensor1.mask),kgrid.x_vec,kgrid.y_vec,kgrid.z_vec);view(127, 18);%disableDefaultInteractivity(ax);

%Plot 2
sensor_xz=zeros(Nx,Nz,kgrid.Nt);
for i=1:kgrid.Nt
    col_i=sensor_data21.p(:,i);
    sensor_i=reshape(col_i,Nx,Nz);
    sensor_xz(:,:,i)=rot90(sensor_i,1);
end
%sensor_xz=(sensor_xz,2);
subplot(1,3,2);flyThrough(sensor_xz,3,100,1,fposxz,(kgrid.x_vec),flipud(kgrid.z_vec),'X-axis','Z-axis');

%Plot 3
subplot(1,3,3);stackedPlot(source3.p((narray*narray+1):end,:));axis square;
xlabel('Time steps');ylabel('Source #');title('Instantaneous Phase')



%--------------Bottom Array Only --------------------
shg;figure('Position', [187,349,1422,558]);  
sgt=sgtitle('Bottom Array Only');sgt.FontSize=20;

%plot 4
subplot(1,3,1);voxelPlot(double(source1.p_mask | sensor1.mask));view(127, 18);

%Plot 5
sensor_xz=zeros(Nx,Nz,kgrid.Nt);
for i=1:kgrid.Nt
    col_i=sensor_data11.p(:,i);
    sensor_i=reshape(col_i,Nx,Nz);
    sensor_xz(:,:,i)=rot90(sensor_i,1);
end
%sensor_xz=(sensor_xz,2);
subplot(1,3,2);flyThrough(sensor_xz,3,50,1,fposxz,(kgrid.x_vec),flipud(kgrid.z_vec),'X-axis','Z-axis');

%Plot 6
subplot(1,3,3);stackedPlot(source3.p(1:(narray*narray),:));axis square;
xlabel('Time steps');ylabel('Source #');title('Instantaneous Phase')


%--------------Both Array --------------------

shg;figure('Position', [187,349,1422,558]);  
sgt=sgtitle('Both Array ON');sgt.FontSize=20;

%plot 7
subplot(1,3,1);voxelPlot(double(source3.p_mask | sensor1.mask));view(127, 18);

%Plot 8
sensor_xz=zeros(Nx,Nz,kgrid.Nt);
for i=1:kgrid.Nt
    col_i=sensor_data31.p(:,i);
    sensor_i=reshape(col_i,Nx,Nz);
    sensor_xz(:,:,i)=rot90(sensor_i,1);
end
%sensor_xz=(sensor_xz,2);
subplot(1,3,2);flyThrough(sensor_xz,3,50,1,fposxz,(kgrid.x_vec),flipud(kgrid.z_vec),'X-axis','Z-axis');

%Plot 9
subplot(1,3,1);voxelPlot(double(source3.p_mask | sensor2.mask));view(127, 18);
sensor_yz=zeros(Ny,Nz,kgrid.Nt);
for i=1:kgrid.Nt
    col_i=sensor_data32.p(:,i);
    sensor_i=reshape(col_i,Ny,Nz);
    sensor_yz(:,:,i)=fliplr(rot90(sensor_i,2));
end
subplot(1,3,3);flyThrough(sensor_yz,3,50,1,fposyz,(kgrid.y_vec),flipud(kgrid.z_vec),'Y-axis','Z-axis');



%--------------RMS Pressure --------------------
shg;figure('Position', [172,372,1438,420]);
sgt=sgtitle('RMS Pressure - Axial and Lateral Sensor');sgt.FontSize=20;
sensor_i=reshape(sensor_data21.p_rms,Nx,Nz);
subplot(1,4,1);imagesc((kgrid.x_vec),flipud(kgrid.z_vec),rot90(sensor_i,1));title('RMS Plot');axis xy;
rectangle('Position',fposxz,'EdgeColor',[0 0 0]);title('Top Array Only - Axial Sensor');
xlabel('X-axis');ylabel('Z-axis');axis square;
lim=clim;c=colorbar;c.Label.String='(Pa)'; c.Label.VerticalAlignment='middle';c.Label.Position=[0.5,-(lim(2)-lim(1))/10];c.Label.Rotation=0;
% c.Label.VerticalAlignment='middle';c.Label.Position=[0.5,0.25];c.Label.HorizontalAlignment='center';
% colormap(cmap);


sensor_i=reshape(sensor_data11.p_rms,Nx,Nz);
subplot(1,4,2);imagesc((kgrid.x_vec),flipud(kgrid.z_vec),rot90(sensor_i,1));title('RMS Plot');axis xy;
rectangle('Position',fposxz,'EdgeColor',[0 0 0]);title('Bottom Array Only - Axial Sensor');
xlabel('X-axis');ylabel('Z-axis');colorbar;axis square;
lim=clim;c=colorbar;c.Label.String='(Pa)'; c.Label.VerticalAlignment='middle';c.Label.Position=[0.5,-(lim(2)-lim(1))/10];c.Label.Rotation=0;
% colormap(cmap);

sensor_i=reshape(sensor_data31.p_rms,Nx,Nz);
subplot(1,4,3);imagesc((kgrid.x_vec),flipud(kgrid.z_vec),rot90(sensor_i,1));title('RMS Plot');axis xy;
rectangle('Position',fposxz,'EdgeColor',[0 0 0]);title('Both Array Combined - Axial Sensor');
xlabel('X-axis');ylabel('Z-axis');colorbar;axis square;
lim=clim;c=colorbar;c.Label.String='(Pa)'; c.Label.VerticalAlignment='middle';c.Label.Position=[0.5,-(lim(2)-lim(1))/10];c.Label.Rotation=0;
% colormap(cmap);

sensor_i=reshape(sensor_data32.p_rms,Nx,Nz);
subplot(1,4,4);imagesc((kgrid.x_vec),flipud(kgrid.z_vec),rot90(sensor_i,1));title('RMS Plot');axis xy;
rectangle('Position',fposyz,'EdgeColor',[0 0 0]);title('Both Array Combined - Lateral Sensor');
xlabel('Y-axis');ylabel('Z-axis');colorbar;axis square;
lim=clim;c=colorbar;c.Label.String='(Pa)'; c.Label.VerticalAlignment='middle';c.Label.Position=[0.5,-(lim(2)-lim(1))/10];c.Label.Rotation=0;
% colormap(cmap);

