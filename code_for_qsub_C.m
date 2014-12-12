% create and run qsub jobs for bubble growth sandbox


if exist('qsubrun_C.out','file')
    system('rm qsubrun_C.out');
end

if exist('qsubrun_C.error','file')
    system('rm qsubrun_C.error');
end


% % first write a one line mfile that calls bubble_growth_sandbox
% fid = fopen(['eval_C_rdot_nuc_rate.m'],'w');
% fprintf(fid,'bubble_growth');
% fclose(fid);

% then write the script to submit the job and call the one liner
fid = fopen(['bubble_growth_C_run.submit'], 'w');
fprintf(fid,['#!/bin/bash\n\n'...
    '#$ -o qsubrun_C.out\n'...
    '#$ -e qsubrun_C.error\n'...
    '#$ -cwd\n'...
    '#$ -S /bin/bash\n'...
    '\nmodule load matlab\n'...
    'matlab -nodesktop < eval_C_rdot_and_nuc_rate.m']);
fclose(fid);


system(['qsub bubble_growth_C_run.submit']);
