function turbulent_res = ReadCklData()

% Reads in data from FINAL_CKL_STATE_Proc# into a convenient Matlab format
% Can also convert to a different number of processors 
n2s = @(str) strrep(num2str(str),'.','');

processed_data_Q=0;
ckl_format = 'float';
base_dir='/Users/jsquire/Documents/MRIDSS/';
data_dir = 'MRIDSS/Data/';
% data_dir = 'ClusterData/';
final_state_dir = 'FINAL_STATE/';

turbulent_res = [];


% for By = [0.05:0.01:0.2]
    By = 0.08;
noise =3;Rm=8000;Pm=2;
run_dir = ['LinearBScans_B' n2s(By)  '_Noise' n2s(noise) 'Rm' n2s(Rm) 'Pm' n2s(Pm) '/']; % file name
% run_dir = ['CheckNoise' n2s(noise) '_Rm' n2s(Rm) 'Pm' n2s(Pm) '/'];
% run_dir = 'CheckNoise5_Rm5000Pm4/';


if ~processed_data_Q
% Some less often changed quantities
    numMFs = 2;
    numFluct = 4;

    P = ReadInputFile([base_dir data_dir run_dir]);
    N = [P.nx P.ny P.nz];
    L = [P.lx P.ly P.lz];

    % Zero processor should be saved, else something wrong
    fid{1} = fopen([base_dir data_dir run_dir final_state_dir 'FINAL_CKL_STATE_Proc0.dat' ]);
    if fid{1} == -1
        error('Unable to open Proc0 state - assuming others do not exist!')
    end

    ii=1;
    while fid{ii}>0  % Loop through finals until one doesn't open to find out number of processes
        ii = ii + 1;
        fid{ii} = fopen([base_dir data_dir run_dir final_state_dir 'FINAL_CKL_STATE_Proc' num2str(ii-1) '.dat' ]);
    end
    % Number of prcoessors given by when it fails
    numprocs = ii-1;

    % Asign cell arrays for the useful things
    Ckl_cell = cell(numprocs,1);
    my_kx_cell = cell(numprocs,1);
    minxy = cell(numprocs,1);
    maxxy = cell(numprocs,1);
    nxy_full = (N(1)-1)*N(2)/2;
    nxy_pproc = nxy_full/numprocs;
    make_complex = @(vec) vec(1,:)+1i*vec(2,:); % Turn fread output complex

    % Read all the data
    for ii=1:numprocs
        mpi_np = fread(fid{ii},1,'int');
        if mpi_np ~= numprocs
            error('Processor numbers dont match up!');
        end
        minxy{ii} = fread(fid{ii},1,'int')+1; % Add one for matlab indexing 
        maxxy{ii} = fread(fid{ii},1,'int')+1;
        % nsteps and t - same on all
        nsteps = fread(fid{ii},1,'int');
        t = fread(fid{ii},1,'double');
        % k_x array - take only the bit belonging to given processor
        my_kx_cell{ii} = fread(fid{ii},[2 nxy_full], 'double');
        my_kx_cell{ii} = make_complex(my_kx_cell{ii});
        my_kx_cell{ii} = my_kx_cell{ii}(minxy{ii}:maxxy{ii});

        % Mean field - same on all 
        MF = fread(fid{ii},[2 numMFs*N(3)], 'double');
        MF =make_complex(MF);
        MF = reshape(MF,[N(3) numMFs]).';

        % Ckl array 
        Ckl_cell{ii} = fread(fid{ii}, [2, nxy_pproc*(numFluct*N(3))^2],ckl_format);
        Ckl_cell{ii} = make_complex(Ckl_cell{ii});
        Ckl_cell{ii} = reshape(Ckl_cell{ii}, [(numFluct*N(3))^2  nxy_pproc]);

        % Check there's no more data
        tmp = fread(fid{ii}, inf);
        if ~isempty(tmp)
            warning('Number of elements did not match up! It is probably all nonsense!!')
        end

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Put Ckl into a more useful form - same as Matlab S3T code
    % This is a little tricky because k arrays get mixed around by remapping

    % Form ky array
    ky = 0:N(2)/2-1;
    ky = repmat(ky,[1, N(1)-1]);
    ky = reshape(ky,[nxy_pproc numprocs]).';
    % Fix kx as read
    read_kx = round(cell2mat(my_kx_cell)/(2i*pi/L(1)));

    % % kx at the current time to match up with read values
    % read_kx = cell2mat(my_kx_cell)/(2i*pi/L(1));
    % kxmat = reshape(kx,[numprocs, nxy_pproc])/(2i*pi/L(1));
    % kymat = reshape(ky,[numprocs, nxy_pproc])/(2i*pi/L(2));
    Cdimz = N(3)*numFluct;
    Ckl = zeros(Cdimz,Cdimz,N(1),N(2)/2);

    nx_array = [0:(N(1)/2-1) -N(1)/2:-1];
    ny_array = 0:(N(2)/2-1);
    rowcol_stor = []; % For error checking
    % Loop through and find matching pair
    for xxx = 1:N(1)
        for yyy = 1:N(2)/2
            % Don't include the nyquist frequency
            if xxx ~= N(1)/2+1
                % Find row (proc number) and column (kx index inside proc)
                [row,col]=find(read_kx==nx_array(xxx) & ky==ny_array(yyy));
                if isempty(row) || isempty(col)
                    warning(['kx = ' num2str(nx_array(xxx)) ', ky = ' num2str(ny_array(yyy)) ...
                        ' not found in read kx data!!']);
                elseif  length(row)>1 || length(col)>1
                    warning(['kx = ' num2str(nx_array(xxx)) ', ky = ' num2str(ny_array(yyy)) ...
                        ' not found more than once!!']);
                else
                    Ckl(:,:,xxx,yyy) = reshape(Ckl_cell{row}(:,col),[Cdimz Cdimz]);
                    rowcol_stor = [rowcol_stor; [row col]]; 
                end

            end
        end
    end

    if length(unique(rowcol_stor,'rows'))~=length(sortrows(rowcol_stor)) || ...
            length(rowcol_stor)~=nxy_full
        warning('Repeated index pair or not enough found!!')
    end

    for ii=1:numprocs
        fclose(fid{ii});
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%    STUDYING THE DATA                   %%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Calculate kx ky kz
    K.KX = (2i*pi/P.lx)*fftshift(-P.nx/2:P.nx/2-1)';
    K.KY = (2i*pi/P.ly)*(0:P.ny/2-1)';
    K.KZ = (2i*pi/P.lz)*fftshift(-P.nz/2:P.nz/2-1)';


    % Energy
    energy=solEnergy(t,P,K,Ckl);
    % Reynolds stress
    [bzuxmuzbx,bzuymuzby] = ReynoldsStresses(t,P,K,Ckl);

    

else
    load([base_dir data_dir run_dir 'ReadCklOutput.mat'])
end



kz1_Bx = reshape(real(bzuxmuzbx(2,:,:)),[P.nx,P.ny/2]);
kz1_By = reshape(real(bzuymuzby(2,:,:)),[P.nx,P.ny/2]);


figure
for kk=1:P.ny/2 % Create legend for each ky
    legendstr{kk} = ['ky = ' num2str(imag(K.KY(kk)))];
end
legendstr{kk+1}='Total';
subplot(211)
plot(fftshift(imag(K.KX)),-fftshift(kz1_Bx,1),fftshift(imag(K.KX)),-sum(fftshift(kz1_Bx,1),2),'k')
title(['Pm = ' num2str(Pm) ' B_y = ' num2str(By) char(10) '  x reynolds stress'])
xlabel('k_x')
legend(legendstr);
subplot(212)
plot(fftshift(imag(K.KX)),fftshift(kz1_By,1),fftshift(imag(K.KX)),sum(fftshift(kz1_By,1),2),'k')
title(['y reynolds stress' char(10) num2str(sum(sum(fftshift(kz1_By))))])
xlabel('k_x')

% turbulent_res = [turbulent_res [By;sum(sum(fftshift(kz1_Bx)));sum(sum(fftshift(kz1_By)))]];

% figure
% [xg,yg,ckl_xy]=xyCorrelationFunction(t,P,K,Ckl);
% contourf(xg,yg,ckl_xy,40,'LineColor','none')


% subplot(313)
% plot(linspace(0,1,length(MF)),real(ifft(MF,[],2)))
% plot(fftshift(imag(K.KX)),fftshift(energy,1))
% xlabel('k_x')
% legend('k_y = 0', ['k_y = ' num2str(imag(K.KY(2)))])
% sum(sum(energy))

% end

end




function P = ReadInputFile(dir_name)
% Reads the input file and returns input parameters
filelist = strsplit(ls(dir_name));
inputfile = [];
for kk=1:length(filelist)
    if ~isempty(strfind(filelist{kk},'.DSSinput'))
        inputfile = filelist{kk};
    end
end
if isempty(inputfile)
    error('Found no input file!');
end
infid = fopen([dir_name inputfile]);
full_file = fscanf(infid,'%c',inf);
full_eqpos = strfind(full_file,'='); % THIS is REALLY INCONVENIENT IN MATLAB



% Get the parameters you want
P.nx = findInputP('nx_','%d');
P.ny = findInputP('ny_','%d');
P.nz = findInputP('nz_','%d');
P.lx = findInputP('lx_','%f');
P.ly = findInputP('ly_','%f');
P.lz = findInputP('lz_','%f');

P.q = findInputP('q_','%f');



fclose(infid);

    function val = findInputP(str, type)
        % Finds the parameter corresponding to str
        vpos = strfind(full_file,str);
        eqpos = full_eqpos(find(full_eqpos>vpos,1));
        val = sscanf(full_file(eqpos+1:end),type,1);
        if isempty(val)
            warning(['Could not find input ' str])
        end
    end

end


function [bzuxmuzbx,bzuymuzby] = ReynoldsStresses(t,P,K,Ckl)
% Calculate the Reynolds stresses for each mode from Ckl to study the
% dominant terms in the dynamos

% To account for having half of a ky grid 
mfac=[1;2*ones(P.ny/2-1,1)];
bzuxmuzbx = zeros(P.nz,P.nx,P.ny/2);
bzuymuzby = zeros(P.nz,P.nx,P.ny/2);

for yyy=1:P.ny/2
    for xxx=1:P.nx
        % Various k values and matrices
        ky=K.KY(yyy);
        kxt=K.KX(xxx)+P.q*t*ky;
        % Form laplacians
        lap2=ky^2 + K.KZ.^2;
        lapF=ky^2+kxt^2+K.KZ.^2;
        if ky==0
            lap2(1)=1;% Shouldn't be energy in this mode anyway
        end
        ilap2 = 1./lap2;
        
        %  Calculate Reynolds stresses as in the matlab code
        ckl=Ckl(:,:,xxx,yyy);

        id=eye(P.nz);
        zs = zeros(P.nz);
        KZ = diag(K.KZ);
        ilap2 = diag(ilap2);
        ubops={[id zs zs zs], [-kxt*ky*ilap2 KZ.*ilap2 zs zs ],...
            [-kxt*KZ.*ilap2 -ky*ilap2 zs zs ],...
            [zs zs id zs], [zs zs -kxt*ky*ilap2 KZ.*ilap2],...
            [zs zs -kxt*KZ.*ilap2 -ky*ilap2]};
        
        % Bx Reynolds stress
        bzuxmuzbx(:,xxx,yyy)=(1/P.nx^2/P.ny^2)*K.KZ.*fft(mfac(yyy)*real(...
            diag(ifft(ifft(  ubops{6}*ckl*(ubops{1}') - ubops{3}*ckl*(ubops{4}')  )')')));
        % By Reynolds stress
        bzuymuzby(:,xxx,yyy)=(1/P.nx^2/P.ny^2)*K.KZ.*fft(mfac(yyy)*real(...
            diag(ifft(ifft(  ubops{6}*ckl*(ubops{2}') - ubops{3}*ckl*(ubops{5}')  )')')));
        
        
    end
end

end


function energy=solEnergy(t,P,K,Ckl)

energy=zeros(P.nx,P.ny/2);
% Calcuates the energy of Ckl solution
if sum(sum(abs(Ckl(:,:,1,1))))>1e-10
    warning('Energy in kx=ky=0 fluctuation tensor')
end

KX = K.KX;
KY=K.KY;
KZ=K.KZ;


mfac=[1;2*ones(P.ny/2-1,1)];

for yyy=1:P.ny/2
    for xxx=1:P.nx
        % Various k values and matrices
        ky=KY(yyy);
        kxt=KX(xxx)+P.q*t*ky;
        % Form laplacians
        lap2=ky^2 + KZ.^2;
        lapF=ky^2+kxt^2+KZ.^2;
        if ky==0
            lap2(1)=1;% Shouldn't be energy in this mode anyway
        end

        Mkl=abs([lapF./lap2;1./lap2;lapF./lap2;1./lap2]);
%         if ts==0
%         disp(['kx = ', num2str(K.kx(xxx)), 'ky = ', num2str(K.ky(yyy))])
%         disp(num2str(mfac(yyy)*sum(Mkl.*diag(Ckl(:,:,xxx,yyy)))));
%         end
        energy(xxx,yyy)=mfac(yyy)*sum(Mkl.*diag(Ckl(:,:,xxx,yyy)));
    end
end
energy=energy/(2*P.nx^2*P.ny^2*P.nz^2);

end


function [xg,yg,ckl_xy]=xyCorrelationFunction(t,P,K,Ckl)
% Returns the xy correlation function as in Simon or Fromang.

% Just uu for now
xyFourier = reshape(Ckl(1,1,:,:),[P.nx,P.ny/2]);% Want the mean component at z1-z2=0
xyFourier = [zeros([P.nx,1]) conj(fliplr(xyFourier(:,2:end))) xyFourier ];

ckl_xy = ifft2(fftshift(fftshift(xyFourier,1)))'; % Transpose so it fits in standard meshgrid format
xg = linspace(-P.lx/2,P.lx/2,P.nx);
yg = linspace(-P.ly/2,P.ly/2,P.ny);
end


