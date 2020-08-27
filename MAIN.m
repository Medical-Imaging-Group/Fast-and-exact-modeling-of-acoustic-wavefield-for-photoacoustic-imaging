
clear all;
rad_r = 22e-3; % sensor radius
num_sensors=60; % number of Sensors
theta = linspace(0, 2*pi, num_sensors+1);
xc = -rad_r.*cos(theta(1:end-1));
yc = rad_r.*sin(-theta(1:end-1));
% Position of the circle center
x0=0;
y0=0;
% x and y coordinates of circle
coords(1:length(xc),1) = xc + x0; 
coords(1:length(xc),2) = yc + y0; 
x_receive=coords(:,2)'; y_receive=coords(:,1)'; % receiver locations

object.Nx = 501;  % number of grid points in the x (row) direction
object.Ny = 501;  % number of grid points in the y (column) direction
object.x = 50.1e-3;              % total grid size [m]
object.y = 50.1e-3;              % total grid size [m]
dx = object.x/object.Nx;          
dy = object.y/object.Ny; 
c_x = ceil(object.Nx/2); c_y = ceil(object.Ny/2);
xx =-100;
yy =-100;
% source locations
x_img = ((1:object.Nx)-((object.Nx+1)/2))*object.x/(object.Nx-1);    x_img = repmat(x_img',1,object.Ny);
y_img = ((1:object.Ny)-((object.Ny+1)/2))*object.y/(object.Ny-1);    y_img = repmat(y_img,object.Nx,1);
M = 100;
N = 100;
indxi = ceil(object.Nx/2) - M  : ceil(object.Nx/2) + M;  
indyi = ceil(object.Ny/2) - N  : ceil(object.Ny/2) + N;  
Nxi = length(indxi);
Nyi = length(indyi);
tl=512; % Considering 512 time samples
Nx = tl*num_sensors;
A_bf = zeros(Nx,Nxi*Nyi);

%Building the system matrix by using all the sources
for i = 1:Nxi
    for j = 1:Nyi
        x = x_img(indyi(j),indxi(i)); y = y_img(indyi(j),indxi(i));
        [tim,mono_fold]=analytic_green([x,y],[x_receive(:),y_receive(:)],1e7,256,1500);
        sd_f=reshape(mono_fold,512,60);
        sd_f=sd_f';
        A_bf(:,(i-1)*Nyi+j) = reshape(sd_f,Nx,1);
    end
end

%%%%%%%%%%%%%%%%%% Blood vessel phantom %%%%%%%%%%%%%%%%%%%
signal_to_noise_ratio = 40;    % [dB]
bv = imread('blood_vessel_grayscale_final.jpg');
BV2=im2bw(bv,0.9);
in1=find(BV2==1);
in0=find(BV2==0);
BV2(in1)=0;
BV2(in0)=1;

%%%%%%%%%%%%%%%%%%%% Derenzo Phantom  %%%%%%%%%%%%%%%%%%
signal_to_noise_ratio = 40;  
bv = imread('derenzo.png');
bv = double(rgb2gray(bv));
ind = find(bv(:)<120);
bv(ind) = 0;
ind = find(bv(:)>120);
bv(ind)=1;
BV1 = imresize(double(bv),[230 230]);
ind = find(BV1(:)<0.9);
BV1(ind) = 0;
ind = find(BV1(:)>0.9);
BV1(ind) = 1;
BV2 = BV1(10:210,20:220);

%%%%%%%%%%%%%%%%%% PAT Phantom %%%%%%%%%%%%%%%%%%%
signal_to_noise_ratio = 40;  
bv = imread('PAT_1.jpg');
bv = double(rgb2gray(bv));
bv = medfilt2(bv,[3 3]);
ind = find(bv(:)<150);
bv(ind) = 0;
ind = find(bv(:)>150);
bv(ind)=1;
BV1 = imresize(double(bv),[210 210],'bicubic');
ind = find(BV1(:)<0.7);
BV1(ind) = 0;
ind = find(BV1(:)>0.7);
BV1(ind) = 1;
BV2 = (BV1(5:205,5:205));            %% 205, 205!! 204 done for BCD
BV2 = medfilt2(BV2,[3 3]);
SE = strel('square',2)
BV2 = imdilate(BV2,SE);
BV3 = edge(BV2,'sobel');
ind = find(BV3(:)>0.5);
BV2(ind)=1;
ind = find(BV2(:)==1);
BV2(ind) = 0.5;
ind = find(BV2(:)==0);
BV2(ind) = 1;
ind = find(BV2(:)==0.5);
BV2(ind) = 0;

medium.sound_speed = 1500;  % [m/s]
time.dt = 5e-8;         % sampling time in sec 5e-8
time.length = 512; 
sensor.mask=[x_receive;y_receive];
object.Nx = 501;  % number of grid points in the x (row) direction
object.Ny = 501;  % number of grid points in the y (column) direction
object.Nx = 501;  % number of grid points in the x (row) direction
object.Ny = 501;  % number of grid points in the y (column) direction
object.x = 50.1e-3;              % total grid size [m]
object.y = 50.1e-3;              % total grid size [m]
dx = object.x/object.Nx;          
dy = object.y/object.Ny; 
object.p0 = zeros(object.Nx, object.Ny);
object.p0(indxi,indyi) = BV2(:,:);
%sensor.frequency_response = [2.25e6 70];
sd2 = forward(object, time, medium, sensor);
sdn2 = addNoise(sd2, signal_to_noise_ratio, 'peak');

%Back projection
x1f= A_bf'*sdn2(:);
x1fbv= 0.04e-6*x1f;
x1fig=reshape(x1fbv,201,201);
figure,imshow(x1fig,[]);colorbar;
