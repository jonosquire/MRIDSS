% MatlabRead.m

n2s = @(str) strrep(num2str(str),'.','');

% Reads in data produced by MRIDSS and plots basic things
base_dir='/Users/jsquire/Documents/MRIDSS/';
data_dir = 'ClusterData/';
data_dir = 'MRIDSS/Data/';
Rm=5000;Pm=1;
By = 0.2;
noise = 50;

run_dir = ['LinearBBxScans_B' n2s(By) 'Rm' n2s(Rm) 'Pm' n2s(Pm) '/']; % file name
run_dir = ['SW_test/']; % file name
Lz = 1;q=1;
%     run_dir = ['LinDyn_uf06bf04/'];


MFdim=32; % Dimension in the mean fields
numMF = 2; % number of mean fields
num_Energy_AM = 2;
num_reynolds = MFdim*2;

fid_en =  fopen([base_dir data_dir run_dir 'energy.dat']);
fid_AM =  fopen([base_dir data_dir run_dir 'angular_momentum.dat']);
% fid_diss =  fopen([base_dir data_dir run_dir 'dissipation.dat']);
fid_rey =  fopen([base_dir data_dir run_dir 'reynolds_stress.dat']);
fid_MF =  fopen([base_dir data_dir run_dir 'mean_fields.dat']);
fid_time = fopen([base_dir data_dir run_dir 'time_vec.dat']);

% Read time vector and simulation length
time = fread(fid_time, inf ,'double');
  time=time(1:end-1);
sim_len= length(time);

figure
% Read in energy and plot
energy = fread(fid_en,[num_Energy_AM,sim_len],'double');
subplot(312);
semilogy(time , sum(energy), 'k',time, energy(:,:) ,'--' );
title('Energy');
if size(energy(1))==4
    legend('Total','Mean U','Mean B','Fluct u','Fluct b')
else
    legend('Total','Mean U','Fluct u')
end

% Angular momentum
AM = fread(fid_AM,[num_Energy_AM,sim_len],'double');
% tAM=AM(1,:)+AM(3,:)-AM(2,:)-AM(4,:);
tAM=AM(1,:)-AM(2,:);
subplot(313);
plot(time,tAM,'k',time, AM,'--' );
title('AM stress');
if size(energy(1))==4
    legend('Total','Mean U','Mean B','Fluct u','Fluct b')
else
    legend('Total','Mean U','Fluct u')
end
% % Dissipation
% diss = fread(fid_diss,[num_Energy_AM,sim_len],'double');
% subplot(312);
% plot(time,sum(diss),'k',time, diss,'--' );
% title('Dissipation');
% legend('Total','Mean U','Mean B','Fluct u','Fluct b')
% 
% % 
% dt=time(2)-time(1);
% disscheck = [-sum(energy);1.5*cumtrapz(time,tAM);-cumtrapz(time,sum(diss))];
% disscheck = [gradient(sum(energy),time);-(1.5*tAM-sum(diss))];
% 
% plot(time,sum(disscheck),'k',time,disscheck)

% %Reynolds stress
% reynolds = fread(fid_rey,[num_reynolds,sim_len],'double');
% endp=length(time);
% plot(time(1:endp), reynolds(1:3,1:endp),time(1:endp),2*reynolds(4:5,1:endp),'LineWidth',1);
% title('Terms in the dynamo equation');
% legend('Shear term','B_y emf','B_y dissipation','10*B_x emf','10*B_x dissipation')
% % plot(time(1:endp), reynolds(1,1:endp),time(1:endp), -50*reynolds(2,1:endp))

% Mean fields;
% tmp = fread(fid_MF , inf ,'double');
% sim_len = length(tmp)/(MFdim*numMF);
tmp = fread(fid_MF , MFdim*numMF*sim_len ,'double');
if length(tmp)~=MFdim*numMF*sim_len
    error('Number of MF elements does not match up');
end
tmp = reshape(tmp,[MFdim*numMF,sim_len]);
z = 0:1/MFdim:(1-1/MFdim);
dz = @(vec) gradient(vec.',z,1).';
dz2 = @(vec) gradient(gradient(vec.',z,1),z,1).';
clear MF;
for ii=0:numMF-1
    MF{ii+1} = tmp(ii*MFdim+1:(ii+1)*MFdim,:);
    MFmax{ii+1} = max(abs(MF{ii+1}));
    MFmean{ii+1} = sqrt(sum(abs(MF{ii+1}).^2)/MFdim);
end
subplot(311);
htime = 1;floor(length(time)/2);
% plot(time(1:htime), log10(MFmean{2}(1:htime)),time(1:htime), log10(MFmean{1}(1:htime)),...
%     time(1:htime), log10(MFmean{3}(1:htime)),time(1:htime), log10(MFmean{4}(1:htime)))
% plot( linspace(0,2*pi,MFdim), MF{2}(:,end),linspace(0,2*pi,MFdim), MF{1}(:,end) )
imagesc(time(htime:end),z ,(MF{2}(:,htime:end)))%,20,'Linestyle','none')
colorbar;
% figure
% title(['Pm = ' num2str(Pm) ' Rm = 8000' num2str(By)])
% for kk=300:10000
%     plot( linspace(0,2*pi,MFdim), MF{2}(:,kk),linspace(0,2*pi,MFdim), MF{1}(:,kk) )
%     title(['t = ' num2str(time(kk))])
% %     ylim([-1e-5,1e-5])
%     drawnow
%     pause(0.02)
% end

fclose(fid_en);
fclose(fid_AM);
fclose(fid_rey);
fclose(fid_MF);
fclose(fid_time);






% checki=find(time>8.99999999,1)
% MFend = [MFend MF{2}(:,checki)];
% diffMF = MFend-repmat(MFend(:,1),1,6);
% diffMF = sqrt(mean(abs(diffMF).^2));
% plot(log10([ 0.002 0.005 0.01  0.02 0.05]),...
%     log10(diffMF(2:end)),'o-')
