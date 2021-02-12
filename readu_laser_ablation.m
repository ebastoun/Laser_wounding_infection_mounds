% readu_laser_ablation - Script to read txt files with displacement vectors strored
% and plot them on top of the Ecadherin images
%
% 1. Read the Ecadherin images and plot on top of them the displacement
% vectors
% 2. Construct kymographs of the average with respect to the wound
% displacement vectors, uy as a function of time, t (min) and vertical position, y (μm)  
% 3. Calculation of mean uy along a distance of 100 μm away from the wound 
% and calculated immediately after ablation
% 4. Storage of the results as a .mat file
% Last modified: Effie Bastounis : 2021-02-10

%First and last frames
kin=2; kfin=107;

%Calibration factor (microns/pixel)
fcal=0.2285;

% location in y (in pixels) where wound starts and ends (assesed based on
% displacement field
win=31; wf=34;

%Directory name where Ecadherin images are stored
dirname='/Users/effiebastounis/Documents/scripts/Pos2/ecad/ecad';

% We used 36x36 windows with overlap of 16 so the resulting matrix with
% displacements is 63x63
sizeu=63;overlap=16;
kymograph=zeros(sizeu,kfin);
umax=10;
outl =5;
s6nl=0;

for k= kin:kfin
    
%Show Ecadherin image    
figure(1)
img = double(imread([dirname sprintf('%3.3d.tif', k)]));
imagesc(log(img));colormap gray;axis image;hold on

%Read and plot on top of the Ecadherin image the displacement vectors
fil1   = '/Users/effiebastounis/Documents/scripts/Pos2/urapiv/ecad';
filename = ([fil1 sprintf('%3.3d.txt', k ) ]);
vec = load(filename); 

n1 = min(find(diff(vec(:,1))<0));
   n2 = length(vec)/n1;

   vec=reshape(vec,n1,n2,5);
   
   x    = vec(:,:,1);
   y    = vec(:,:,2);
   u    = vec(:,:,3);
   v    = vec(:,:,4);
   s6n = vec(:,:,5);

   if s6nl>0

      % signal to noise check
      
      s6n_low = find(s6n<s6nl);
      
      w = 3; rad = 2;
      kernel = zeros(2*w+1);
      for i=1:2*w+1;
          kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/rad^2);
      end
      kernel(w+1,w+1)=0;
      kernel = kernel/sum(sum(kernel));
   
      u(s6n_low) = 0;
      v(s6n_low) = 0;
      
      tmpv = (conv2(v,kernel,'same'));
      tmpu = (conv2(u,kernel,'same'));
      
      % Let's throw the outlayers out:

      u(s6n_low) = tmpu(s6n_low); 
      v(s6n_low) = tmpv(s6n_low); 
      u(s6n<s6nl)=NaN;
      v(s6n<s6nl)=NaN;

     uu=sqrt(u.^2+v.^2);
     u(uu>umax)=NaN;
     v(uu>umax)=NaN;
  end 
   
  if outl>0

      % Adaptive Local Median filtering
   
      w = 2; rad = 1;
      kernel = zeros(2*w+1);
      for i=1:2*w+1;
          kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/rad^2);
      end
      kernel(w+1,w+1)=0;
      kernel = kernel/sum(sum(kernel));
   
      tmpv = (conv2(v,kernel,'same'));
      tmpu = (conv2(u,kernel,'same'));
   
      lmtv_p = mean(mean(tmpv(2:end-1,2:end-1))) + ...
             outl*std(reshape(tmpv(2:end-1,2:end-1),(n1-2)*(n2-2),1));
      lmtv_m = mean(mean(tmpv(2:end-1,2:end-1))) - ...
             outl*std(reshape(tmpv(2:end-1,2:end-1),(n1-2)*(n2-2),1));
      lmtu_p = mean(mean(u(2:end-1,2:end-1))) + ...
             outl*std(reshape(u(2:end-1,2:end-1),(n1-2)*(n2-2),1));
      lmtu_m = mean(mean(u(2:end-1,2:end-1))) - ...
             outl*std(reshape(u(2:end-1,2:end-1),(n1-2)*(n2-2),1));

      u_out_p = find(u>lmtu_p);
      u_out_m = find(u<lmtu_m);
      v_out_p = find(v>lmtv_p);
      v_out_m = find(v<lmtv_m);
   
      % Let's throw the outlayers out:

      u(u_out_m) = tmpu(u_out_m); 
      u(v_out_m) = tmpu(v_out_m); 
      v(u_out_m) = tmpv(u_out_m); 
      v(v_out_m) = tmpv(v_out_m); 

      u(u_out_p) = tmpu(u_out_p); 
      u(v_out_p) = tmpu(v_out_p); 
      v(u_out_p) = tmpv(u_out_p); 
      v(v_out_p) = tmpv(v_out_p); 

  end  

%Increase value of qq if you want to plot less vectors   
qq=1;

%Increase/decrease value if you want to make the vectors look bigger/smaller   
size=50;
hold on
quiver(x(1:qq:end,1:qq:end),y(1:qq:end,1:qq:end),u(1:qq:end,1:qq:end)*size*fcal,v(1:qq:end,1:qq:end)*size*fcal,'AutoScale','off','color','g');  
axis off;set(gcf,'inverthardcopy','off')    
title(sprintf('Frame %d',k));set(gca,'Fontsize',14); axis image;%pause;

% Save Ecahderin image with displacement vectors supeimposed
eval(['print -dtiffnocompression /Users/effiebastounis/Documents/scripts/Pos2/urapiv/frame-arrows' int2str(k)]);

% Calculate real values of u and v, plot and save meam displacement along the y
% axis (vertical to the wound)
ur=(u)'*fcal;vr=(v)'*fcal;
figure;plot([1:sizeu]*overlap*fcal,mean(vr,2),'b');title(sprintf('Frame %d',k));
xlabel('Distance along y-axis (um)');ylabel('Mean displacement (um)');
set(gca,'FontSize',18,'DefaultAxesFontName', 'Arial');
eval(['print -dtiffnocompression /Users/effiebastounis/Documents/scripts/Pos2/urapiv/integral' int2str(k)]);

% Plot and save 2D map of v (uy), displacement along the direction vertical to the wound
figure;imagesc([1:sizeu]*overlap*fcal,[1:sizeu]*overlap*fcal,vr);title(sprintf('Frame %d, uy (um)',k));colorbar
set(gca,'FontSize',18,'DefaultAxesFontName', 'Arial');
eval(['print -dtiffnocompression /Users/effiebastounis/Documents/scripts/Pos2/urapiv/vrmap' int2str(k)]);

% Construct the kymograph of mean uy (v)
kymom=mean(vr,2);
kymographm(:,k)=kymom; close all;
end 

% Plot final kymograph excluding noisy vectors in the wound area
figure
kymographfm=kymographm;

% Get rid of the noisy vectors in the wound area
kymographfm(wf:end,:)=-kymographfm(wf:end,:);

% Change sign of displacement values below the wound so that displacements
% towards the wound are positive (red) and away from the wound negative
% (blue)
kymographfm(win:wf,:)=0;
imagesc([1:kfin],[1:sizeu]*overlap*fcal,kymographfm);
set(gca,'FontSize',18,'DefaultAxesFontName', 'Arial');
xlabel('Time post-wounding (min)');ylabel('Distance (um)');colorbar
caxis([-0.1 0.1]); colormap redblue
eval(['print -dtiffnocompression /Users/effiebastounis/Documents/scripts/Pos2/urapiv/kymographm' int2str(k)]);

% Calculate mean recoil displacement (rd) taking into account 100 um away from the wound
% Given 100 um/(fcal*overlap)= 27 pixels
% You can change also average (or not) the first 3 frames post-wounding as below
rd_mounders= mean(mean(kymographfm(wf+1:wf+27,kin:kin+2)));
rd_surrounders= mean(mean(kymographfm(win-28:win-1, kin:kin+2)));

save /Users/effiebastounis/Documents/scripts/Pos2/urapiv/deformations.mat
