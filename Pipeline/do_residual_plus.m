function [rDir,rMag,u_mean,v_mean,w_mean]=do_residual_plus(u,v,w,dt)
% DO_RESIDUAL Takes the u and v vectors of a model output and calculates
% the long-term direction and magnitude for that data.
% 
%   [RDIR,RMAG,URES,VRES]=DO_RESIDUAL(U,V,DT) takes the residual direction
%   (RDIR) and magnitude RMAG) of the data in U and V sampled at interval
%   DT. URES and UDIR are the summed U and V positions (the raw data for a
%   progresive vector diagram). Direction output is in degrees, vector
%   magnitude in units/s.
% 
% Pierre Cazenave PML 20/03/2012.
% 

% Loosely based on my original dfsuResidual.m and processResidual function
% for DHI's MIKE21 software, which in turn were based on Dave Lambkin's
% residual analysis scripts.
% 
% TODO: Make it possible to specify the average for all layers (i.e. NZ is
% all layers). 

% Let's do it...

%%i think dt has to be in days, so what is the time interval that the
%%layers are at, so in my case it's 1/24 because there are 24 layers, each
%%an hour? or should it be 1/(24*60*60) so that it's in seconds?? if my
%%input is all dx/dt, dy/dt etc, then dt is 1/(24*60*60)
%looking at the dfsuResidual.m file, he used number of timesteps and then
%timestep in seconds. Here, nTimesteps is 24 because it's the layers i have
%so timestep would be 60*60 if it's in seconds?
%  duration=timestep*num_timesteps;

toSecFactor=24*60*60;

nElements=size(u,1);
nLayers=size(u,2);
nTimeSteps=size(u,3);

% Some tidal assumptions. This will need to change in areas in which the
% diurnal tide dominates over the semidiurnal. 
tideCycle=(12+(25/60))/24; %units is days? 0.51 days/tide cycle
tideWindow=ceil(tideCycle/dt); %this should be in hours because it is used to index the layers later? so now, 13 if dt=1/24
%tideDuration=(mean((dt*nTimeSteps)-tideCycle)-mean(tideCycle))*toSecFactor; %turning into seconds from days
tideDuration=(dt*nTimeSteps)*toSecFactor;


% Preallocate outputs.
uRes=zeros(nElements,nLayers,nTimeSteps);
vRes=zeros(nElements,nLayers,nTimeSteps);
wRes=zeros(nElements,nLayers,nTimeSteps);
uSum=nan(nElements,nTimeSteps,nLayers); %not sure why timesteps are in the middle now...
vSum=nan(nElements,nTimeSteps,nLayers);
wSum=nan(nElements,nTimeSteps,nLayers);
% uStart=nan(nElements,nLayers);
% vStart=nan(nElements,nLayers);
% wStart=nan(nElements,nLayers);
% uEnd=nan(nElements,nLayers);
% vEnd=nan(nElements,nLayers);
% wEnd=nan(nElements,nLayers);
u_mean=nan(nElements,nLayers);
v_mean=nan(nElements,nLayers);
w_mean=nan(nElements,nLayers);



for hh=1:nLayers %for each layer
    uSum(:,:,hh)=cumsum(squeeze(u(:,hh,:)),2); %%squeeze gives array where rows are velocity and columns are time, so you have all data for each depth
    vSum(:,:,hh)=cumsum(squeeze(v(:,hh,:)),2); %%cumsum(,2) takes cumulative sum along rows, so result is array of nelementsxtimesteps, giving cumulative sum at one point over time
    wSum(:,:,hh)=cumsum(squeeze(w(:,hh,:)),2); %% wSum's first page now holds the result of cumulative sum at the first depth over time
    for ii=1:nTimeSteps %for each timestep
        uRes(:,hh,ii)=uRes(:,hh,ii)+(uSum(:,ii,hh));%*(dt*toSecFactor)); %%multiplying here by 24*60*60 which puts it into hours???
        vRes(:,hh,ii)=vRes(:,hh,ii)+(vSum(:,ii,hh));%.*(dt*toSecFactor)); %% .* is element-wise multiplication
        wRes(:,hh,ii)=wRes(:,hh,ii)+(wSum(:,ii,hh));%.*(dt*toSecFactor)); %%uRes's pages are timesteps, so it gives 
        %so now uRes etc has cumulative velocity in m/hour 
    end
    % uStart(:,hh)=mean(squeeze(uRes(:,hh,1:tideWindow)),2); %so tideWindow has to be a whole number!
    % vStart(:,hh)=mean(squeeze(vRes(:,hh,1:tideWindow)),2);
    % wStart(:,hh)=mean(squeeze(wRes(:,hh,1:tideWindow)),2);
    % uEnd(:,hh)=mean(squeeze(uRes(:,hh,end-tideWindow:end)),2);
    % vEnd(:,hh)=mean(squeeze(vRes(:,hh,end-tideWindow:end)),2);
    % wEnd(:,hh)=mean(squeeze(wRes(:,hh,end-tideWindow:end)),2);
end

u_mean=uRes(:,:,nTimeSteps)/nTimeSteps;%/tideDuration;
v_mean=vRes(:,:,nTimeSteps)/nTimeSteps;%/tideDuration;
w_mean=wRes(:,:,nTimeSteps)/nTimeSteps;%/tideDuration;

% uDiff=uEnd-uStart;
% vDiff=vEnd-vStart;
% wDiff=wEnd-wStart; %this is net movement over the whole tide duration

% new_u=uDiff/tideDuration;
% new_v=vDiff/tideDuration;
% new_w=wDiff/tideDuration;

% Calculate direction and magnitude.
% rDir=atan2(uDiff,vDiff)*(180/pi); % in degrees.
% rMag=sqrt(uDiff.^2+vDiff.^2+wDiff.^2)/tideDuration; % in units/s.
rDir=atan2(v_mean,u_mean)*(180/pi); % in degrees.
rMag=sqrt(u_mean.^2+v_mean.^2+w_mean.^2); % in units/s.
%%dividing by number of seconds to turn into units/s
%which means it was m/tide duration

end