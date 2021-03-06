function WriteMRIDSSinput

% Writes an input file for the c++ MRIDSS program - useful to run a series
% with different parameters

% Optionally, MRIDSScan also be run from here

n2s = @(str) strrep(num2str(str),'.','');

MRIdir='../MRIDSS/';% Base directory

for Bval = [ 1.0]; 
    for B = {'Bx'}
    
        
    Rm=5000;
    Pm=1;
        filename =  ['Linear' B{1} 'Scan_B' n2s(Bval) 'Rm' n2s(Rm) 'Pm' n2s(Pm)]; % file name;
    % Simulation type
    equation_type = 'MHD_fullBQlin';
    % Parameters, comment out some of these for different parameter scans
    
    L=[1,2,1];
    N=[32 32 48];
    % Physical parameters
    q=1.5;
    omega = q*2/3;
    % Mean vertical field
    B0z = 0.0;

    nu = Pm/Rm; % viscosity
    eta = 1/Rm; % resistivity
    f_noise = 5; % Driving noise
    noise_range = [0 20000];
    % Time steps
    dt = 0.01; % Timestep
    CFL = 1.7; % CFL number
    time_interval = [0 30]; % Time interval
    tv_save = 0.2; % Save energy, momentum, MFs etc. every tv_save
    full_save = time_interval(2)+1;
    % Flags for saving
    save_enQ=1; save_amQ=1; save_dissQ=1;
    save_ReyQ = 1; %Reynolds stress
    save_mfQ=1; % Mean-fields
    % Initial conditions
    if strcmp(B{1},'Bx');bxoby = -0.1;else bxoby=1;end
    init_By = bxoby*Bval;
    % Simulation flags
    remapQ =1; % Remapping
    QuasiLinearQ = 0; % Reynolds stress feedback
    
    StartFromSavedQ=0;



    % Check for input files already there
    file_list = ls(MRIdir);
    if ~isempty(strfind(file_list,'.DSSinput'))
        warning('MRIDSS directory already contains an input folder: moving file to .OldInputFiles!');
        file_list = strsplit(file_list);
        for kk = 1:length(file_list)
            if ~isempty(strfind(file_list{kk},'.DSSinput'))
                movefile([MRIdir file_list{kk}],[MRIdir '.OldInputFiles/']);
            end
        end
    end


    fid = fopen([MRIdir filename '.DSSinput'],'w+t');
    if fid<0
        error('Failed to open file!!');
    end

    fprintf(fid,'Input file for MRIDSS - Generated automatically by WriteMRIDSSinput.m\n\n','char');
    fprintf(fid,'Format: variable_ = var (can have arbitrary number of spaces around = )\n\n', 'char');
    
    fprintf(fid,'// Simulation type - used as a check\n','char');
    fprintf(fid,['equations_to_use_ = ' equation_type ' \n' ],'char');
    
    % Printing variables
    fprintf(fid, '// Box parameters\n','char');
    fprintf(fid, 'lx_ = %15.15f, ly_ = %15.15f, lz_ = %15.15f \n',L(1),L(2),L(3));
    fprintf(fid, 'nx_ = %d,  ny_ = %d,  nz_ = %d \n\n',N(1),N(2),N(3));

    fprintf(fid, '// Time parameters\n','char');
    fprintf(fid, 'dt_ = %15.15f\n',dt);
    fprintf(fid, 'CFL_ = %15.2f\n',CFL);
    fprintf(fid, 't_initial_ = %15.15f,   t_final_ = %15.15f,\n\n',time_interval(1),time_interval(2));

    fprintf(fid, '// Saving\n','char');
    fprintf(fid, 'tv_save_ = %15.15f,   // Save energy, momentum etc. time-step\n',tv_save);
    fprintf(fid, 'full_save_ = %15.15f,   // Save full solution\n',full_save);
    fprintf(fid, 'save_energy? =%d  save_angular_mom? =%d     save_dissipation?=%d, // On-off flags\n',save_enQ,save_amQ,save_dissQ);
    fprintf(fid, 'save_reynolds_stress?= %d \n', save_ReyQ);
    fprintf(fid, 'mean_field_save? = %d  // Save mean fields\n\n',save_mfQ);

    fprintf(fid, '// General physical parameters\n','char');
    fprintf(fid, 'q_ = %15.15f \n',q);
    fprintf(fid, 'omega_ = %15.15f \n',omega);
    fprintf(fid, 'B0z_ = %15.15f \n',B0z);
    fprintf(fid, 'nu_ = %15.15f,  eta_ = %15.15f,   // Dissipation parameters \n', nu,eta);
    fprintf(fid, 'f_noise_ = %15.15f \n\n' , f_noise);
    fprintf(fid, 'noise_range_low_ = %15.15f,  noise_range_high_ = %15.15f\n\n' , noise_range(1),noise_range(2));
    
    fprintf(fid, '// Initial conditions\n','char');
    fprintf(fid, 'initial_By_ = %15.15f \n',init_By);

    fprintf(fid, '// Flags for simulation\n','char');
    fprintf(fid, 'Remap? = %d \n',remapQ);
    fprintf(fid, 'QuasiLinear? = %d \n\n',QuasiLinearQ);
    
    fprintf(fid,'// Restart simulation\n','char');
    fprintf(fid,'start_from_saved? = %d \n' ,StartFromSavedQ);


    fclose(fid);

    disp(['Done writing ' filename '.DSSinput'])
    disp(' ')
    if isempty(strfind(getenv('PATH'),'/usr/local/bin'))
        setenv('PATH',[getenv('PATH') ':/usr/local/bin']);
    end
    !mpiexec -np 8 ~/Documents/MRIDSS/mridss_prog

    disp(' ')
    disp(['Done running case ' filename ])
    disp(' ')
    disp('<<<<<<<<<<<<<<>>>>>>>>>>>>>')
    disp('<<<<<<<<<<<<<<>>>>>>>>>>>>>')
    disp(' ')
    disp(' ')
    end
    

end






