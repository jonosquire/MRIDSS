% MatlabRead_extractGrowth.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Loop through a series of runs and extract growth rate or saturation


for Pmkk = [1 2 3 4];
    tmp_array = [];
for Rm = 1000:500:5000;

% Reads in data produced by MRIDSS and plots basic things
base_dir='/Users/jsquire/Documents/MRIDSS/';
data_dir = 'ClusterData/RmScan/';
run_dir = ['RmScanNoise3_Rm' num2str(Rm) 'Pm' num2str(Pmkk) '/'];
MFdim=48; % Dimension in the mean fields
numMF = 2; % number of mean fields
num_Energy_AM = 4;

fid_en =  fopen([base_dir data_dir run_dir 'energy.dat']);
fid_AM =  fopen([base_dir data_dir run_dir 'angular_momentum.dat']);
fid_rey =  fopen([base_dir data_dir run_dir 'reynolds_stress.dat']);
fid_MF =  fopen([base_dir data_dir run_dir 'mean_fields.dat']);
fid_time = fopen([base_dir data_dir run_dir 'time_vec.dat']);

% Read time vector and simulation length
time = fread(fid_time, inf ,'double');
sim_len= length(time);

% Mean fields
tmp = fread(fid_MF , inf ,'double');
if length(tmp)~=MFdim*numMF*sim_len
    error('Number of NF elements does not match up');
end
tmp = reshape(tmp,[MFdim*numMF,sim_len]);
for ii=0:numMF-1
    MF{ii+1} = tmp(ii*MFdim+1:(ii+1)*MFdim,:);
    MFmax{ii+1} = max(abs(MF{ii+1}));
end

% % Angular momentum
% energy = fread(fid_en,[num_Energy_AM,sim_len],'double');
% plot(time , sum(energy), 'k',time, energy ,'--' );
% 
% title(run_dir);
% satval = ginput(1);
% satval=satval(2)
% tmp_array = [tmp_array [noisekk;  satval]];


htime = floor(length(time)/2);
plot(time(1:htime), log10(MFmax{2}(1:htime)))
title(run_dir);

timebnds = ginput(2);
tb =[find(time>timebnds(1,1),1) find(time>timebnds(2,1),1)]; sort(tb);
growth_rate = (log(MFmax{2}(tb(2)))-log(MFmax{2}(tb(1))))/(time(tb(2))-time(tb(1)));

tmp_array = [tmp_array [noisekk;  growth_rate]];

% satval = ginput(1);
% satval=10^(satval(2))
% tmp_array = [tmp_array [noisekk;  satval]];


fclose(fid_en);
fclose(fid_AM);
fclose(fid_rey);
fclose(fid_MF);
fclose(fid_time);

end
GrowthRate.(['Pm' strrep(num2str(Pmkk),'.','')]) = tmp_array;
end