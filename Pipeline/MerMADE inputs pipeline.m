%This script takes hydrodynamic files, interpolates them to user-defined
%resolution and breaks into user-defined number of layers for MerMADE
%inputs

%%
%Set up the script

%Tell the script where your files are. Matlab should have moved to the
%folder where you opened the script from, but just in case:
base_dir= 'C:/Users/Rey/Documents/Scallops/'; %!!! change this filepath to where you have your fvcom toolbox saved
cd(base_dir); %this changes the working directory to that filepath
%add the "libraries" needed
addpath([base_dir 'fvcom-toolbox/fvcom_prepro'])
addpath([base_dir 'fvcom-toolbox/fvcom_postproc'])
addpath([base_dir 'fvcom-toolbox/utilities'])
addpath([base_dir 'RPSstuff-master/RPSstuff-master/RPSstuff']) %for the ll2utm funcion 
addpath([base_dir 'm_map1.4/m_map/']) %for the m_proj function within ll2utm

%define path and filename to work with
daily_path= 'Data/';
month = 'April/'; %!!! change this to whichever month and day the data is coming from, or whatever you've named that folder

%note that this month's .nc file should be the only one in the folder
findfile=strcat([base_dir daily_path month],'*.nc');
dinfo = dir(findfile);
filename=fullfile(dir(findfile).folder, dir(findfile).name);

%if you want to have a look at what's in the .nc file, remove the % in the
%next line
%ncdisp(filename);

%setting your parameters

%!!! latitude and longitude of your map
%the script will subset the .nc file data to these specifications
min_long=-7.1;
max_long=-4.8;
min_lat=53.9;
max_lat=55.4;

%resolution info (in metres)
xy_res=500;
z_res=10;
max_depth=160;
nlayers=max_depth/z_res; 
%when layers are created, it counts from 0 so even though this is 16, there will be 17 layers created

%%
%Read in the mesh structure and time array
M.time=ncread(filename, 'time');
M.nv=ncread(filename,'nv');
%lat/lon for the elements
M.lon= ncread(filename, 'lon');
M.lat = ncread(filename, 'lat');
M.lonc = ncread(filename, 'lonc'); %degrees east
M.latc= ncread(filename, 'latc'); %degrees north
M.sf_depth=ncread(filename, 'h_center'); % seafloor depth below geoid
M.siglay=ncread(filename, 'siglay');
M.siglay_center=ncread(filename, 'siglay_center');

%hydro data in m/s
M.u = ncread(filename, 'u'); %these are calculated based on elements NOT nodes
M.v = ncread(filename, 'v');
M.w = ncread(filename, 'ww');

%calculating residuals
[res_dir, res_spd, new_u,new_v,new_w]=do_residual_plus(M.u,M.v,M.w,1/24);
M.res_dir=res_dir;
M.res_spd=res_spd*60*60; %turning it into m/hr
M.u_res=new_u;
M.v_res=new_v;
M.w_res=new_w;


%FVCOM reference time is Modified Julian time 1858 -11 -17 00:00:00
%Convert time to matlab time
M.time_matlab= double(M.time) + datetime(1858,11,17);
M.lonc(M.lonc>180)=M.lonc(M.lonc>180)-360;
M.lon(M.lon>180)=M.lon(M.lon>180)-360;

%if you want to Examine the current direction and speed in a plot, remove
%the % symbols from the following lines until line 91
% Quiver vector (plot every other depth average velocity vector on unstructured grid):
% Q.X = M.lonc(1:15:end);
% Q.Y = M.latc(1:15:end);
% Q.U = M.u_res(1:15:end,1);
% Q.V = M.v_res(1:15:end,1);
% Q.scale= 4;
% Q.colour = "white";
% plot_fvcom_field(M, M.res_spd(:,1), 'pll','qui', Q);
% axis([min_long max_long min_lat max_lat])


%%
indx=0;
for r=1:size(M.lonc)
    if M.lonc(r) >= min_long
        if M.lonc(r) <= max_long
         if M.latc(r)  >= min_lat 
             if M.latc(r) <=max_lat
                 if indx==0
                     indx=r;
                 else
                     indx(end+1)=r;
                 end
           
             end
         end
        end
    end
    
end

for i=1:length(indx)
    fprintf('indx[%d] \n',i);
    
     if i==1
        S.lonc_sub=M.lonc(indx(i));
        S.latc_sub=M.latc(indx(i));
        S.nv=M.nv(indx(i));
        S.u=M.u_res(indx(i), :);
        S.v=M.v_res(indx(i), :);
        S.w=M.w_res(indx(i), :);
        S.depth=M.sf_depth(indx(i));
        S.siglay_c=M.siglay_center(indx(i),:);
        
     else
        S.lonc_sub(end+1)=M.lonc(indx(i));
        S.latc_sub(end+1)=M.latc(indx(i));
        S.nv(end+1,:)=M.nv(indx(i));
        S.u(end+1,:,1)=M.u_res(indx(i), :);
        S.v(end+1, :,1)=M.v_res(indx(i), :);
        S.w(end+1,:,1)=M.w_res(indx(i), :);
        S.depth(end+1)=M.sf_depth(indx(i));
        S.siglay_c(end+1,:)=M.siglay_center(indx(i),:);
     end
end


save([daily_path month 'cropped_extent'], "S");
%!!! if you want to come back to this point without starting from scratch,
%remove the % in the next line and read in the file created up to this
%point
%S=load([daily_path month '\cropped_extent.mat'], "-mat", "S").S;

%scatteredInterpolant
X1= double(S.lonc_sub(:));
Y1= double(S.latc_sub(:));
[xlo,yla] = ll2utm(X1, Y1, 30);
S.lonc_m=xlo;
S.latc_m=yla;

% calculate absolute depth in each layer
for e=1:length(S.depth) %go through each element
    for l=1:width(S.siglay_c) %for each layer
%         fprintf('doing element %d , layer %d \n', e,l);
        S.Seafloor_depth(e,l)=-1*S.siglay_c(e,l)*S.depth(e);
    end
    fprintf('doing element %d \n', e);
end


save([daily_path '\' month 'cropped_extent'], "S");
%!!! if you want to come back to this point without starting from scratch,
%remove the % in the next line and read in the file created up to this
%point
%S= load([daily_path month '\cropped_extent.mat'], "-mat", "S").S;


%now i have my cropped data, i need to interpolate for each variable
%the new struct with the interpolated data will have:
% lonc_m, latc_m, nv, u, v, w, depth

%I need to create a matrix for x, y, and z for scattered interpolant, where
%each is latc_m*how many layers the dataset has i want long and 1 column (this is a column vector)
x_col= (repelem(S.lonc_m(:,1),width(S.siglay_c))); %doing this per layer, so i'll have each lonc_m value repeated x times
y_col=(repelem(S.latc_m(:,1),width(S.siglay_c)));

for i=1:length(S.latc_m)
    fprintf("doing element %d \n", i);
    if i==1
        z_col=S.Seafloor_depth(i,:)';
        
    else
        z_col=[z_col;S.Seafloor_depth(i,:)'];
    end

end

%now I need a U, V, W value that is also latc_m*nlay long
for i=1:length(S.latc_m)
    fprintf("doing element %d \n", i);
    if i==1
        u_int=S.u(i,:)';
        v_int=S.v(i,:)';
        w_int=S.w(i,:)';
        
    else
        u_int=[u_int;S.u(i,:)'];
        v_int=[v_int;S.v(i,:)'];
        w_int=[w_int;S.w(i,:)'];
    end

end

u_data=[x_col, y_col, z_col, u_int];
v_data=[x_col, y_col, z_col, v_int];
w_data=[x_col, y_col, z_col, w_int];
all_data=[x_col, y_col,z_col, u_int, v_int, w_int];
save([daily_path '\' month 'all_matrix'], "all_data");
save([daily_path '\' month 'u_matrix'], "u_data");
save([daily_path '\' month 'v_matrix'], "v_data");
save([daily_path '\' month 'w_matrix'], "w_data");

%!!! if you want to come back to this point without starting from scratch,
%remove the % in the next line and read in the file created up to this
%point
%u_data=load([daily_path month 'u_matrix.mat']).u_data;
%v_data=load([daily_path month 'v_matrix.mat']).v_data;
%w_data=load([daily_path month 'w_matrix.mat']).w_data;
%all_data=load([daily_path month 'all_matrix.mat']).all_data;

x_col=double(all_data(:,1));
y_col=double(all_data(:,2));
z_col=double(all_data(:,3));
u_int=double(all_data(:,4));
v_int=double(all_data(:,5));
w_int=double(all_data(:,6));

 fprintf("doing interpolation of u");
[x,y,z]=meshgrid(min(S.lonc_m):xy_res:max(S.lonc_m), max(S.latc_m):-1*xy_res:min(S.latc_m), 0:z_res:max_depth);
si_u=scatteredInterpolant(double(x_col), double(y_col), double(z_col),u_int, 'nearest');
fint_u= si_u(x,y,z);
save([daily_path '\' month 'u_fint'], "fint_u");
 fprintf("doing interpolation of v");
si_v=scatteredInterpolant(double(x_col), double(y_col), double(z_col),v_int, 'nearest');
fint_v=si_v(x,y,z);
save([daily_path '\' month 'v_fint'], "fint_v");
 fprintf("doing interpolation of w");
si_w=scatteredInterpolant(double(x_col), double(y_col), double(z_col),w_int, 'nearest');
fint_w=si_w(x,y,z);
save([daily_path '\' month 'w_fint'], "fint_w");
 fprintf("doing interpolation of depth");
%For seafloor depth, i only really need 1 layer
si_depth=scatteredInterpolant(double(S.lonc_m(:,1)), double(S.latc_m(:,1)), double(S.depth(1,:)'), "nearest");
fint_sfdepth= si_depth(x,y);
save([daily_path '\' month 'sfdepth_fint'], "fint_sfdepth");

%creating an object of just the interpolant data
R.x_m=min(S.lonc_m):xy_res:max(S.lonc_m);
R.y_m=max(S.latc_m):-1*xy_res:min(S.latc_m);
R.z_m=0:z_res:max_depth;
R.u_fint= fint_u;
R.v_fint= fint_v;
R.w_fint= fint_w;
R.sf_depth=fint_sfdepth;

save([daily_path '\' month 'study_area_' string(xy_res)], "R");
%!!! if you want to come back to this point without starting from scratch,
%remove the % in the next line and read in the file created up to this
%point
%R=load([daily_path month 'study_area_' string(xy_res)]).R;

%create layers
 fprintf("creating layers");
for la=0:nlayers
    u_lay= R.u_fint(:,:,la+1);
    f_name= sprintf('%s_%d.txt',"u",la);
    writematrix(u_lay, [daily_path '\' month f_name], "Delimiter", "tab");
    
    v_lay= R.v_fint(:,:,la+1);
    f_name= sprintf('%s_%d.txt',"v",la);
    writematrix(v_lay, [daily_path '\' month f_name], "Delimiter", "tab");
    
    w_lay= R.w_fint(:,:,la+1);
    f_name= sprintf('%s_%d.txt',"w",la);
    writematrix(w_lay, [daily_path '\' month f_name], "Delimiter", "tab"); %this still writes as scientific notation because values are so small
    
end

sfd_lay=R.sf_depth(:,:,1);
f_name= 'sf_depth.txt';
writematrix(sfd_lay, [daily_path '\' month f_name], "Delimiter", "tab");

