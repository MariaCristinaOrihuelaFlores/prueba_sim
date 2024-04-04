clear; clc; close all;
addpath("ustb");
addpath("k-Wave");
%%
DATA_CAST = 'gpuArray-single'; %gpuArray- 
RUN_SIMULATION  = true;
%%
f0 = 8e6;       % pulse center frequency [Hz]
cycles = 4;       % number of cycles in pulse
c0 = 1540;      % medium speed of sound [m/s]
rho0 = 1000;    % medium density [kg/m3]
N = 1;            % number of plane waves 
prb = uff.linear_array();
prb.N = 32;                  % number of elements
prb.pitch = 300e-6;           % probe pitch in azimuth [m]
prb.element_width = 300e-6;   % element width [m]
prb.element_height = 5000e-6; %y-axis
%fig_handle = prb.plot([],'Linear array');
%% 
dx=prb.pitch/120;
PML_size = 100; % size of the PML in grid points
Nx=round(10e-3/dx); Nx=Nx+mod(Nx,2);
Nz=round(30e-3/dx); Nz=Nz+mod(Nz,2); 
grid_width=Nx*dx;
grid_depth=Nz*dx;
domain=uff.linear_scan('x_axis', linspace(-grid_width/2,grid_width/2,Nx).', 'z_axis', linspace(0,grid_depth,Nz).');
kgrid = kWaveGrid(domain.N_z_axis, domain.z_step, domain.N_x_axis, domain.x_step);
%% 
Nz_vacio = round(20e-3/dx);
vacio_region = zeros(Nz_vacio,Nx);
rad_inclusion = 2e-3;
small_diameter = 150e-6;
large_diameter = 150e-6;
lat_pos = grid_width/2;
ax_pos = 25e-3;
num_small = 1000;%calculate_num_circ(radius_inclusion,small_concentration,large_diameter);
num_large = 100;%calculate_num_background(grid_depth-1e-3,grid_width,radius_inclusion,large_concentration,large_diameter);
%%
scattering_region = define_map_region(Nz-Nz_vacio, Nx, num_large, round(large_diameter/(2*dx)), round(rad_inclusion/dx), [round(ax_pos/dx)-8001,round(lat_pos/dx)]);
region_all = vertcat(vacio_region,scattering_region);
%% 
%figure;imagesc(domain.x_axis*1e3,domain.z_axis*1e3,region_all);axis equal tight;
%%
%savepath=[pwd,'\sim\'];
%save(fullfile(savepath,'med_attemp4.mat'),'region_all'); 
%% 
%load med_attemp4.mat

scattering_region = region_all;
scattering_region_copy = scattering_region;

scattering_region(find(scattering_region == 0)) = 1540;
scattering_region(find(scattering_region == 1)) = 1950;
medium.sound_speed = scattering_region;

% density
scattering_region_copy(find(scattering_region_copy == 0)) = 1000;
scattering_region_copy(find(scattering_region_copy == 1)) = 960;
medium.density = scattering_region_copy;
%% 
% figure;
% subplot(1,2,1);
% imagesc(domain.x_axis*1e3,domain.z_axis*1e3,medium.sound_speed); colormap gray; colorbar; axis equal tight;
% xlabel('x [mm]');
% ylabel('z [mm]');
% title('c_0 [m/s]');
% subplot(1,2,2);
% imagesc(domain.x_axis*1e3,domain.z_axis*1e3,medium.density); colormap gray; colorbar; axis equal tight;
% xlabel('x [mm]');
% ylabel('z [mm]');
% title('\rho [kg/m^3]');
%% 
cfl = 0.3; % depth
t_end = 2*sqrt(grid_width.^2+grid_depth.^2)/mean(medium.sound_speed(:));
kgrid.makeTime(medium.sound_speed,cfl,t_end);
%%
angles = 0;
seq =uff.wave();
seq.apodization = uff.apodization('f_number',1,'window',uff.window.hamming,'focus',uff.scan('xyz',[0 0 10e-3]));
seq.source.distance = -Inf;
seq.probe = prb;
seq.sound_speed = 1540;    % reference speed of sound [m/s]
seq.delay = min(seq.delay_values);
seq.source.azimuth = angles;
%%
source_pixels={};
element_sensor_index = {};
n=1;
% element_position: prb.x
for m=1:prb.N_elements
%             start points    last points        
    %plot((prb.x(m)+[-prb.width(m)/2 prb.width(m)/2])*1e3,[0 0],'k+-'); hold on; grid on;
    %pause(2)
    source_pixels{m}=find(abs(domain.x-prb.x(m))<prb.width(m)/2 & abs(domain.y-prb.y(m))<prb.height(m) & abs(domain.z-prb.z(m))<=domain.z_step/2);
    element_sensor_index{m} = n:n+numel(source_pixels{m})-1;
    n=n+numel(source_pixels{m});
end
% sensor mask
sensor.mask = zeros(domain.N_z_axis, domain.N_x_axis);
for m=1:prb.N_elements
    sensor.mask(source_pixels{m}) = sensor.mask(source_pixels{m}) + 1;
end

% source mask
source.u_mask=sensor.mask;
%%
delay=seq.delay_values-seq.delay;
denay=round(delay/kgrid.dt);
seq.delay = seq.delay - cycles/f0/2;
    
% offsets
tone_burst_offset = [];
for m=1:prb.N_elements
    tone_burst_offset = [tone_burst_offset repmat(denay(m),1,numel(source_pixels{m}))];
end
source.ux = toneBurst(1/kgrid.dt, f0, cycles, 'SignalOffset', tone_burst_offset);   % create the tone burst signals
source.uy = 0.*source.ux;
source.u_mode ='dirichlet';
    
% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size,'DataCast',DATA_CAST, 'PlotPML', false, 'Smooth', false};
    
% run the simulation
sensor_data(:,:,1) = permute(kspaceFirstOrder2DG(kgrid, medium, source, sensor, input_args{:}),[2 1]);
%% 
sensor_data(isnan(sensor_data))=0;
element_data=zeros(numel(kgrid.t_array),prb.N_elements,numel(seq));
for m=1:prb.N_elements
    if  ~isempty(element_sensor_index{m})
        element_data(:,m,:)=bsxfun(@times,sqrt(1./kgrid.t_array).',trapz(kgrid.y(source_pixels{m}),sensor_data(:,element_sensor_index{m},:),2));
    end
end
channel_data = uff.channel_data();
channel_data.probe = prb;
channel_data.sequence = seq;
channel_data.initial_time = 0;
channel_data.sampling_frequency = 1/kgrid.dt;
channel_data.data = element_data;
channel_data.data(isnan(channel_data.data))=0;

savepath=[pwd,'\channel\'];
save(fullfile(savepath,'prueba_attemp.mat'),'channel_data');