% create and run qsub jobs for bubble growth sandbox


if exist('qsubrun.out','file')
    system('rm qsubrun.out');
end

if exist('qsubrun.error','file')
    system('rm qsubrun.error');
end


% first write a one line mfile that calls bubble_growth_sandbox
fid = fopen(['bubble_growth_input.m'],'w');
fprintf(fid,'bubble_growth');
fclose(fid);

% then write the script to submit the job and call the one liner
fid = fopen(['bubble_growth_run.submit'], 'w');
fprintf(fid,['#!/bin/bash\n\n'...
    '#$ -o qsubrun.out\n'...
    '#$ -e qsubrun.error\n'...
    '#$ -cwd\n'...
    '#$ -S /bin/bash\n'...
    '\nmodule load matlab\n'...
    'matlab -nodesktop < bubble_growth_input.m']);
fclose(fid);


system(['qsub bubble_growth_run.submit']);
