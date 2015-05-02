% fiddle with constants

m_start = 0;

% cd('sweeping_constants')


C_death_rate = logspace(-3,6,20);
C_erf_death_rate = logspace(-2,2,20);
C_coalescence = logspace(-3,6,20);

% C_rdot = 1;
% C_nuc_rate = 1;
% C_death_rate = 1;
% C_erf_death_rate = 1;


m = 0;
for j = 1:length(C_coalescence)
    for k = 1:length(C_death_rate)
        for n = 1:length(C_erf_death_rate)
            m = m + 1;
            
            if m > m_start
                
                [stat, n_l_s] = system('qstat | grep -c ^');
                n_lines = str2double(n_l_s);
                
                while n_lines > 470
                    pause(5*60)
                    [stat, n_l_s] = system('qstat | grep -c ^');
                    n_lines = str2double(n_l_s);
                end
                
                
                %                     error_filename = ['qsubrun' num2str(m) '.error'];
                %                     % check if the error file exists
                %                     if exist(error_filename,'file')
                %
                %                         % if it does, see if it's size = 0
                %                         % if it does, then it's already been run with no error
                %                         s = dir(error_filename);
                %                         if s.bytes == 0
                %                             % no error when it was already run
                %                             run_this_one = 0;
                %                         else
                %                             % some error occured. run again
                %                             run_this_one = 1;
                %                         end
                %                     else
                %                         % error file doesn't exist -> run it
                run_this_one = 1;
                %                     end
                
                
                
                if run_this_one
                    
                    save_filename = ['bubble_sim_data_number_' num2str(m) '.mat'];
                    
                    system(['./code_for_qsub_separate_jobs.sh ' num2str(m) ' ' ...
                        num2str(C_coalescence(j)) ' ' ...
                        num2str(C_death_rate(k)) ' ' ...
                        num2str(C_erf_death_rate(n)) ' ' save_filename]);
                    %             pause
                    
                end
            end
            
        end
        
        
    end
end

