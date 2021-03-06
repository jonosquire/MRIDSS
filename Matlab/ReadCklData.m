function  ReadCklData()

% Reads in data from FINAL_CKL_STATE_Proc# into a convenient Matlab format
% Can also convert to a different number of processors 
n2s = @(str) strrep(num2str(str),'.','');

processed_data_Q=0;
ckl_format = 'float';
base_dir='/Users/jsquire/Documents/MRIDSS/';
data_dir = 'MRIDSS/Data/';
% data_dir = 'Clus,terData/RmScan/';
final_state_dir = 'FINAL_STATE/';

turbulent_res = [];

sp=0;
store_full_rey = struct;
store_total_rey = struct;
for Pm = [1]
    for Rm = [ 1]
        storetmp = [];
        for By=0.05

% run_dir = ['LinearBA0521Scans_B' n2s(By) 'Rm' n2s(Rm) 'Pm' n2s(Pm) '/']; % file name
        run_dir =  ['Linear' 'Bx' 'Scan_B' n2s(0.1) 'Rm' n2s(5000) 'Pm' n2s(1) '/']; % file name;


% 

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
        MF = make_complex(MF);
        MF = reshape(MF,[N(3) numMFs]).';

        % Ckl array 
        Ckl_cell{ii} = fread(fid{ii}, [2, nxy_pproc*(numFluct*N(3))^2],ckl_format);
        Ckl_cell{ii} = make_complex(Ckl_cell{ii});
        try 
          Ckl_cell{ii} = reshape(Ckl_cell{ii}, [(numFluct*N(3))^2  nxy_pproc]);
        catch
            Ckl_cell{ii} = zeros([(numFluct*N(3))^2  nxy_pproc]);
            warning(['Did not work at Pm' n2s(Pm) 'Rm' n2s(Rm) 'By' n2s(By) '\n nz = ' n2s(N(3))])
        end

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
    read_kx_init = round(cell2mat(my_kx_cell)/(2i*pi/L(1)));
    read_kx = mod(read_kx_init+N(1)/2-1,N(1)-1)-(N(1)/2-1); % With new remapping kx gets smaller continuously, 
    % this mod shouldn't affect read_kx from old data
    

    % % kx at the current time to match up with read values
    % read_kx = cell2mat(my_kx_cell)/(2i*pi/L(1));
    % kxmat = reshape(kx,[numprocs, nxy_pproc])/(2i*pi/L(1));
    % kymat = reshape(ky,[numprocs, nxy_pproc])/(2i*pi/L(2));
    Cdimz = N(3)*numFluct;
    Ckl = zeros(Cdimz,Cdimz,N(1)-1,N(2)/2);

    nx_array = [0:(N(1)/2-1) -N(1)/2+1:-1];
    ny_array = 0:(N(2)/2-1);
    kx_array = zeros(N(1)-1,N(2)/2); % Store the original version of kx as read from file
    rowcol_stor = []; % For error checking
    % Loop through and find matching pair
    for xxx = 1:N(1)-1
        for yyy = 1:N(2)/2
            % Find row (proc number) and column (kx index inside proc)
            [row,col]=find(read_kx==nx_array(xxx) & ky==ny_array(yyy));
            if isempty(row) || isempty(col)
                warning(['kx = ' num2str(nx_array(xxx)) ', ky = ' num2str(ny_array(yyy)) ...
                    ' not found in read kx data!!']);
            elseif  length(row)>1 || length(col)>1
                warning(['kx = ' num2str(nx_array(xxx)) ', ky = ' num2str(ny_array(yyy)) ...
                    ' found more than once!!']);
            else
                Ckl(:,:,xxx,yyy) = reshape(Ckl_cell{row}(:,col),[Cdimz Cdimz]);
                kx_array(xxx,yyy) = read_kx_init(row,col);
                rowcol_stor = [rowcol_stor; [row col]]; 
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
    K.KY = (2i*pi/P.ly)*(0:P.ny/2-1)';
    K.kx_array = bsxfun(@plus,(2i*pi/L(1))*kx_array,P.q*t*K.KY.');
    K.KZ = (2i*pi/P.lz)*fftshift(-P.nz/2:P.nz/2-1)';

    % Put the ks in order
    [K.kx_array, reord]=sort(imag(K.kx_array),1);
    K.kx_array = 1i*K.kx_array;
    for kk=1:P.ny/2
        Ckl(:,:,:,kk)=Ckl(:,:,reord(:,kk)',kk);
    end

    % Energy
     [energy,energyfull]=solEnergy(P,K,Ckl);
    % Reynolds stress
    [bzuxmuzbx,bzuymuzby] = ReynoldsStresses(t,P,K,Ckl);

    

else
    load([base_dir data_dir run_dir 'ReadCklOutput.mat'])
end



kz1_Bx = reshape(real(K.KZ(2,:,:).*bzuxmuzbx(2,:,:)),[P.nx-1,P.ny/2]);
kz1_By = reshape(real(K.KZ(2,:,:).*bzuymuzby(2,:,:)),[P.nx-1,P.ny/2]);

% kyarr=repmat(imag(K.KY)',[P.nx-1,1]);
% kz1_Bx(kyarr==0)=0;
% kz1_By(kyarr==0)=0;
% energy(kyarr==0)=0;


figure
cmap = colormap;cvec = 0:1/64:1-1/64;
plotnum = P.ny/2;
for kk=1:plotnum % Create legend for each ky
    legendstr{kk} = ['ky = ' num2str(imag(K.KY(kk)))];
    colorlist{kk} = interp1(cvec,cmap,(kk-1)/plotnum);
end
legendstr{kk+1}='Total';

sp=sp+1;
totB = [sum(sum(kz1_Bx)) ;sum(sum(kz1_By));sum(sum(energy))];
for yyy=1:plotnum
    subplot(3,1,sp)
    hold on
    plot(imag(K.kx_array(:,yyy)),kz1_Bx(:,yyy),'o-','Color',colorlist{yyy},'MarkerEdgeColor',colorlist{yyy},'MarkerFaceColor',colorlist{yyy});
    title(['Pm = ' num2str(Pm) ' B_y = ' num2str(By) '  x reynolds stress:  '  num2str(sum(sum(kz1_Bx)))])
    xlabel('k_x')
    subplot(312)
    hold on
    plot(imag(K.kx_array(:,yyy)),kz1_By(:,yyy),'o-','Color',colorlist{yyy},'MarkerEdgeColor',colorlist{yyy},'MarkerFaceColor',colorlist{yyy});
    title(['y reynolds stress:  '  num2str(sum(sum(kz1_By)))])
    xlabel('k_x')
    subplot(313)
    hold on
    plot(imag(K.kx_array(:,yyy)),energy(:,yyy),'o-','Color',colorlist{yyy},'MarkerEdgeColor',colorlist{yyy},'MarkerFaceColor',colorlist{yyy});
    title(['Energy'])
    xlabel('k_x')
end
legend(legendstr)
sum(sum(energy))

% storetmp = [storetmp totB];
% store_full_rey.(['Pm' n2s(Pm) 'Rm' n2s(Rm) 'By' n2s(By)]).Bx = kz1_Bx;
% store_full_rey.(['Pm' n2s(Pm) 'Rm' n2s(Rm) 'By' n2s(By)]).Bx = kz1_By;
% store_full_rey.(['Pm' n2s(Pm) 'Rm' n2s(Rm) 'By' n2s(By)]).K = K;
% disp(['Done: Pm' n2s(Pm) 'Rm' n2s(Rm) 'By' n2s(By)]);

        end % By loop 
        store_total_rey.(['Pm' n2s(Pm) 'Rm' n2s(Rm)]) = storetmp;
    end % Rm loop 
end % Pm loop
% store_total_rey.By = 0.02:0.02:0.2;
% 
% R = store_total_rey;
% save([base_dir data_dir 'FullData_wBxA0521'],'store_full_rey')
% save([base_dir data_dir 'DataTotals_wBxA0521'],'R')


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
bzuxmuzbx = zeros(P.nz,P.nx-1,P.ny/2);
bzuymuzby = zeros(P.nz,P.nx-1,P.ny/2);

kxfac = (2i*pi/P.lx);
for yyy=1:P.ny/2
    for xxx=1:P.nx-1
        % Various k values and matrices
        ky=K.KY(yyy);
%         kxt=K.KX(xxx)+P.q*t*ky;
%         kxt = kxfac*(mod(kxt/kxfac+P.nx/2,P.nx)-P.nx/2);
        kxt = K.kx_array(xxx,yyy);
        % Fix for the new remapping style, mod back to -nx/2->nx/2
        
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
        bzuxmuzbx(:,xxx,yyy)=(1/P.nx^2/P.ny^2).*fft(mfac(yyy)*real(...
            diag(ifft(ifft(  ubops{6}*ckl*(ubops{1}') - ubops{3}*ckl*(ubops{4}')  )')')));
        % By Reynolds stress
        bzuymuzby(:,xxx,yyy)=(1/P.nx^2/P.ny^2).*fft(mfac(yyy)*real(...
            diag(ifft(ifft(  ubops{6}*ckl*(ubops{2}') - ubops{3}*ckl*(ubops{5}')  )')')));
        
        
    end
end

end


function [energy,energyfull]=solEnergy(P,K,Ckl)

energy=zeros(P.nx-1,P.ny/2);
% Calcuates the energy of Ckl solution

mfac=[1;2*ones(P.ny/2-1,1)];

for yyy=1:P.ny/2
    for xxx=1:P.nx-1
        % Various k values and matrices
        ky=K.KY(yyy);
        kxt = K.kx_array(xxx,yyy);
        % Form laplacians
        lap2=ky^2 + K.KZ.^2;
        lapF=ky^2+kxt^2+K.KZ.^2;
        if ky==0
            lap2(1)=1;% Shouldn't be energy in this mode anyway
        end
        on = ones(size(lap2));on(1)=1;

        Mkl=abs([on.*lapF./lap2;on.*1./lap2;0*lapF./lap2;0*1./lap2]);
        
%         if ts==0
%         disp(['kx = ', num2str(K.kx(xxx)), 'ky = ', num2str(K.ky(yyy))])
%         disp(num2str(mfac(yyy)*sum(Mkl.*diag(Ckl(:,:,xxx,yyy)))));
%         end
        energyfull(:,xxx,yyy)=Mkl.*diag(Ckl(:,:,xxx,yyy));
        energy(xxx,yyy)=mfac(yyy)*sum(Mkl.*diag(Ckl(:,:,xxx,yyy)));
        
        if kxt==0.0 && ky==0.0
            if sum(sum(abs(Ckl(:,:,xxx,yyy))))>1e-10
                warning('Energy in kx=ky=0 fluctuation tensor')
            end
        end

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




