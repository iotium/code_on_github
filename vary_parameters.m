% parameters to vary
% number of nodes (spatial)
% number of moments
% integration tolerance

% qdot_lw
% coalescence
% nucleation rate
% nucleation distribution


% first, number of nodes:

N_nodes = [8 16 32 64 128 256 512];
N_mom = [4 6 8 10];
rel_tol = 10.^[-1 -2 -3 -4 -5 -6];
N_rw = [8 16 32 64 128 256];

C_qdot_lw = [ 0.125 0.25 0.5 2 4 8]*2e-5;
C_coalescence = [0.125 0.25 0.5 2 4 8]*0.1;
% C_coalescence = logspace(-5,0,6);
C_coalescence = [C_coalescence(:) ones(length(C_coalescence),2)];
C_nuc_rate = [0.125 0.25 0.5 2 4 8]*1e-8;

inputs.rel_tol = 1e-3;
inputs.N_mom = 4;
inputs.N_nodes = 32;
inputs.N_rw = 32;

inputs.C_qdot_lw = 2e-5;
inputs.C_coalescence = [0.1 1 1]; % collision efficiency, laminar shear, turbulence
inputs.C_nuc_rate = 1e-8;


% ------ vary N_rw
k = 200;
for i = 1:length(N_rw)
    % running qsub jobs
    
    m_file_name = ['qsub_mfile_' num2str(k) '.m'];
    fid = fopen(m_file_name,'w');
    
    fprintf(fid, ['inputs.N_nodes = %0.4g\n'], inputs.N_nodes);
    fprintf(fid, ['inputs.N_mom = %0.4g\n'], inputs.N_mom);
    fprintf(fid, ['inputs.rel_tol = %0.4g\n'], inputs.rel_tol);
    fprintf(fid, ['inputs.N_rw = %0.4g\n'], N_rw(i));

    fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], inputs.C_qdot_lw);
    fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], inputs.C_coalescence);
    fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], inputs.C_nuc_rate);
    fprintf(fid, ['bubble_growth_1D(inputs)']);
    
    fclose(fid);
    
    submit_file_name = ['submit_file_' num2str(k) '.submit'];
    output_file_name = ['output_file_' num2str(k) '.out'];
    error_file_name = ['error_file_' num2str(k) '.error'];
    
    fid = fopen(submit_file_name, 'w');
    
    fprintf(fid,['#!/bin/bash\n\n']);
    fprintf(fid,['#$ -o %s\n'], output_file_name);
    fprintf(fid,['#$ -e %s\n'], error_file_name);
    fprintf(fid,['#$ -cwd\n'...
        '#$ -S /bin/bash\n'...
        '\nmodule load matlab\n']);
    fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
    fclose(fid);
    
    system(['qsub ' submit_file_name ]);
    
    k = k + 1;
end

% % ------ vary C_qdot_lw
% 
% k = 100;
% for i = 1:length(C_qdot_lw)
%     % running qsub jobs
%     
%     m_file_name = ['qsub_mfile_' num2str(k) '.m'];
%     fid = fopen(m_file_name,'w');
%     
%     fprintf(fid, ['inputs.N_nodes = %0.4g\n'], inputs.N_nodes);
%     fprintf(fid, ['inputs.N_mom = %0.4g\n'], inputs.N_mom);
%     fprintf(fid, ['inputs.rel_tol = %0.4g\n'], inputs.rel_tol);
%         fprintf(fid, ['inputs.N_rw = %0.4g\n'], inputs.N_rw);
%     fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], C_qdot_lw(i));
%     fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], inputs.C_coalescence);
%     fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], inputs.C_nuc_rate);
%     fprintf(fid, ['bubble_growth_1D(inputs)']);
%     
%     fclose(fid);
%     
%     submit_file_name = ['submit_file_' num2str(k) '.submit'];
%     output_file_name = ['output_file_' num2str(k) '.out'];
%     error_file_name = ['error_file_' num2str(k) '.error'];
%     
%     fid = fopen(submit_file_name, 'w');
%     
%     fprintf(fid,['#!/bin/bash\n\n']);
%     fprintf(fid,['#$ -o %s\n'], output_file_name);
%     fprintf(fid,['#$ -e %s\n'], error_file_name);
%     fprintf(fid,['#$ -cwd\n'...
%         '#$ -S /bin/bash\n'...
%         '\nmodule load matlab\n']);
%     fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
%     fclose(fid);
%     
%     system(['qsub ' submit_file_name ]);
%     
%     k = k + 1;
% end
% 
% % ------ vary C_coalescence
% 
% for i = 1:length(C_qdot_lw)
%     % running qsub jobs
%     
%     m_file_name = ['qsub_mfile_' num2str(k) '.m'];
%     fid = fopen(m_file_name,'w');
%     
%     fprintf(fid, ['inputs.N_nodes = %0.4g\n'], inputs.N_nodes);
%     fprintf(fid, ['inputs.N_mom = %0.4g\n'], inputs.N_mom);
%     fprintf(fid, ['inputs.rel_tol = %0.4g\n'], inputs.rel_tol);
%         fprintf(fid, ['inputs.N_rw = %0.4g\n'], inputs.N_rw);
% 
%     fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], inputs.C_qdot_lw);
%     fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], C_coalescence(i,:));
%     fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], inputs.C_nuc_rate);
%     fprintf(fid, ['bubble_growth_1D(inputs)']);
%     
%     fclose(fid);
%     
%     submit_file_name = ['submit_file_' num2str(k) '.submit'];
%     output_file_name = ['output_file_' num2str(k) '.out'];
%     error_file_name = ['error_file_' num2str(k) '.error'];
%     
%     fid = fopen(submit_file_name, 'w');
%     
%     fprintf(fid,['#!/bin/bash\n\n']);
%     fprintf(fid,['#$ -o %s\n'], output_file_name);
%     fprintf(fid,['#$ -e %s\n'], error_file_name);
%     fprintf(fid,['#$ -cwd\n'...
%         '#$ -S /bin/bash\n'...
%         '\nmodule load matlab\n']);
%     fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
%     fclose(fid);
%     
%     system(['qsub ' submit_file_name ]);
%     
%     k = k + 1;
% end
% 
% % ------ vary C_nuc_rate
% 
% for i = 1:length(C_nuc_rate)
%     % running qsub jobs
%     
%     m_file_name = ['qsub_mfile_' num2str(k) '.m'];
%     fid = fopen(m_file_name,'w');
%     
%     fprintf(fid, ['inputs.N_nodes = %0.4g\n'], inputs.N_nodes);
%     fprintf(fid, ['inputs.N_mom = %0.4g\n'], inputs.N_mom);
%     fprintf(fid, ['inputs.rel_tol = %0.4g\n'], inputs.rel_tol);
%         fprintf(fid, ['inputs.N_rw = %0.4g\n'], inputs.N_rw);
% 
%     fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], inputs.C_qdot_lw);
%     fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], inputs.C_coalescence);
%     fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], C_nuc_rate(i));
%     fprintf(fid, ['bubble_growth_1D(inputs)']);
%     
%     fclose(fid);
%     
%     submit_file_name = ['submit_file_' num2str(k) '.submit'];
%     output_file_name = ['output_file_' num2str(k) '.out'];
%     error_file_name = ['error_file_' num2str(k) '.error'];
%     
%     fid = fopen(submit_file_name, 'w');
%     
%     fprintf(fid,['#!/bin/bash\n\n']);
%     fprintf(fid,['#$ -o %s\n'], output_file_name);
%     fprintf(fid,['#$ -e %s\n'], error_file_name);
%     fprintf(fid,['#$ -cwd\n'...
%         '#$ -S /bin/bash\n'...
%         '\nmodule load matlab\n']);
%     fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
%     fclose(fid);
%     
%     system(['qsub ' submit_file_name ]);
%     
%     k = k + 1;
% end


% % ------ vary N_nodes
% 
% k = 1;
% for i = 1:length(N_nodes)
%     % running qsub jobs
%     
%     m_file_name = ['qsub_mfile_' num2str(k) '.m'];
%     fid = fopen(m_file_name,'w');
%     
%     fprintf(fid, ['inputs.N_nodes = %0.4g\n'], N_nodes(i));
%     fprintf(fid, ['inputs.N_mom = %0.4g\n'], inputs.N_mom);
%     fprintf(fid, ['inputs.rel_tol = %0.4g\n'], inputs.rel_tol);
%     fprintf(fid, ['inputs.N_rw = %0.4g\n'], inputs.N_rw);

%     fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], inputs.C_qdot_lw);
%     fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], inputs.C_coalescence);
%     fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], inputs.C_nuc_rate);
%     fprintf(fid, ['bubble_growth_1D(inputs)']);
%     
%     fclose(fid);
%     
%     submit_file_name = ['submit_file_' num2str(k) '.submit'];
%     output_file_name = ['output_file_' num2str(k) '.out'];
%     error_file_name = ['error_file_' num2str(k) '.error'];
%     
%     fid = fopen(submit_file_name, 'w');
%     
%     fprintf(fid,['#!/bin/bash\n\n']);
%     fprintf(fid,['#$ -o %s\n'], output_file_name);
%     fprintf(fid,['#$ -e %s\n'], error_file_name);
%     fprintf(fid,['#$ -cwd\n'...
%         '#$ -S /bin/bash\n'...
%         '\nmodule load matlab\n']);
%     fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
%     fclose(fid);
%     
%     system(['qsub ' submit_file_name ]);
%     
%     k = k + 1;
% end
% 
% 
% % ------ vary N_mom
% 
% for i = 1:length(N_mom)
%     % running qsub jobs
%     
%     m_file_name = ['qsub_mfile_' num2str(k) '.m'];
%     fid = fopen(m_file_name,'w');
%     
%     fprintf(fid, ['inputs.N_nodes = %0.4g\n'], inputs.N_nodes);
%     fprintf(fid, ['inputs.N_mom = %0.4g\n'], N_mom(i));
%     fprintf(fid, ['inputs.rel_tol = %0.4g\n'], inputs.rel_tol);
%     fprintf(fid, ['inputs.N_rw = %0.4g\n'], inputs.N_rw);

%     fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], inputs.C_qdot_lw);
%     fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], inputs.C_coalescence);
%     fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], inputs.C_nuc_rate);
%     fprintf(fid, ['bubble_growth_1D(inputs)']);
%     
%     fclose(fid);
%     
%     submit_file_name = ['submit_file_' num2str(k) '.submit'];
%     output_file_name = ['output_file_' num2str(k) '.out'];
%     error_file_name = ['error_file_' num2str(k) '.error'];
%     
%     fid = fopen(submit_file_name, 'w');
%     
%     fprintf(fid,['#!/bin/bash\n\n']);
%     fprintf(fid,['#$ -o %s\n'], output_file_name);
%     fprintf(fid,['#$ -e %s\n'], error_file_name);
%     fprintf(fid,['#$ -cwd\n'...
%         '#$ -S /bin/bash\n'...
%         '\nmodule load matlab\n']);
%     fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
%     fclose(fid);
%     
%     system(['qsub ' submit_file_name ]);
%     
%     k = k + 1;
% end
% 
% % ------ vary rel_tol
% 
% for i = 1:length(rel_tol)
%     % running qsub jobs
%     
%     m_file_name = ['qsub_mfile_' num2str(k) '.m'];
%     fid = fopen(m_file_name,'w');
%     
%     fprintf(fid, ['inputs.N_nodes = %0.4g\n'], inputs.N_nodes);
%     fprintf(fid, ['inputs.N_mom = %0.4g\n'], inputs.N_mom);
%     fprintf(fid, ['inputs.rel_tol = %0.4g\n'], rel_tol(i));
%     fprintf(fid, ['inputs.N_rw = %0.4g\n'], inputs.N_rw);

%     fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], inputs.C_qdot_lw);
%     fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], inputs.C_coalescence);
%     fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], inputs.C_nuc_rate);
%     fprintf(fid, ['bubble_growth_1D(inputs)']);
%     
%     fclose(fid);
%     
%     submit_file_name = ['submit_file_' num2str(k) '.submit'];
%     output_file_name = ['output_file_' num2str(k) '.out'];
%     error_file_name = ['error_file_' num2str(k) '.error'];
%     
%     fid = fopen(submit_file_name, 'w');
%     
%     fprintf(fid,['#!/bin/bash\n\n']);
%     fprintf(fid,['#$ -o %s\n'], output_file_name);
%     fprintf(fid,['#$ -e %s\n'], error_file_name);
%     fprintf(fid,['#$ -cwd\n'...
%         '#$ -S /bin/bash\n'...
%         '\nmodule load matlab\n']);
%     fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
%     fclose(fid);
%     
%     system(['qsub ' submit_file_name ]);
%     
%     k = k + 1;
% end
% 
% % ------ vary rel_tol
% 
% for i = 1:length(rel_tol)
%     % running qsub jobs
%     
%     m_file_name = ['qsub_mfile_' num2str(k) '.m'];
%     fid = fopen(m_file_name,'w');
%     
%     fprintf(fid, ['inputs.N_nodes = %0.4g\n'], inputs.N_nodes);
%     fprintf(fid, ['inputs.N_mom = %0.4g\n'], inputs.N_mom);
%     fprintf(fid, ['inputs.rel_tol = %0.4g\n'], rel_tol(i));
%     fprintf(fid, ['inputs.N_rw = %0.4g\n'], inputs.N_rw);

%     fprintf(fid, ['inputs.C_qdot_lw = %0.4g\n'], inputs.C_qdot_lw);
%     fprintf(fid, ['inputs.C_coalescence = [%0.4g %0.4g %0.4g]\n'], inputs.C_coalescence);
%     fprintf(fid, ['inputs.C_nuc_rate = %0.4g\n'], inputs.C_nuc_rate);
%     fprintf(fid, ['bubble_growth_1D(inputs)']);
%     
%     fclose(fid);
%     
%     submit_file_name = ['submit_file_' num2str(k) '.submit'];
%     output_file_name = ['output_file_' num2str(k) '.out'];
%     error_file_name = ['error_file_' num2str(k) '.error'];
%     
%     fid = fopen(submit_file_name, 'w');
%     
%     fprintf(fid,['#!/bin/bash\n\n']);
%     fprintf(fid,['#$ -o %s\n'], output_file_name);
%     fprintf(fid,['#$ -e %s\n'], error_file_name);
%     fprintf(fid,['#$ -cwd\n'...
%         '#$ -S /bin/bash\n'...
%         '\nmodule load matlab\n']);
%     fprintf(fid,['matlab -nodesktop < %s'], m_file_name);
%     fclose(fid);
%     
%     system(['qsub ' submit_file_name ]);
%     
%     k = k + 1;
% end

% write a short m file that sets inputs and calls bubble growth

% write a .submit file that calls this file



% for i = 1:length(C_nuc_rate);
%     inputs.C_nuc_rate = C_nuc_rate(i);
%     try
%         bubble_growth_1D(inputs)
%     end
% end

% for i = 1:length(C_coalescence(:,1));
%     inputs.C_coalescence = C_coalescence(i,:);
%     try
%         bubble_growth_1D(inputs)
%     end
% end

% for i = 1:length(C_qdot_lw);
%     inputs.C_qdot_lw = C_qdot_lw(i);
%     try
%         bubble_growth_1D(inputs)
%     end
% end


% bubble_growth_1D(inputs);

%

%
% inputs.rel_tol = 1e-4;
% inputs.N_nodes = 50;
%
% for i = 1:length(N_mom)
%     inputs.N_mom = N_mom(i);
%     try
%     bubble_growth_1D(inputs);
%     end
% end
%
% inputs.N_nodes = 50;
% inputs.N_mom = 6;
%
% for i = 1:length(rel_tol)
%     inputs.rel_tol = rel_tol(i);
%     try
%     bubble_growth_1D(inputs);
%     end
% end