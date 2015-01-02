% fiddle with constants


C_rdot = linspace(1,10,7);
C_nuc_rate = logspace(0,6,15);
C_death_rate = logspace(0,6,15);
C_erf_death_rate = logspace(-2,1,8);

% C_rdot = 1;
% C_nuc_rate = 1;
% C_death_rate = 1;
% C_erf_death_rate = 1;


m = 0;
n_submitted = 0;
for i = 1:length(C_rdot)
    for j = 1:length(C_nuc_rate)
        for k = 1:length(C_death_rate)
            for n = 1:length(C_erf_death_rate)
                m = m + 1;
                
                %     constants.C_rdot = varargin{1};
                %     constants.C_nuc_rate = varargin{2};
                %     constants.n_nuc_freq = 3;
                %     constants.checking_gauss_error = 0;
                %     constants.C_death_rate = varargin{3};
                %     constants.C_erf_death_rate = varargin{4};
                %     constants.coalescence_switch = 'on';
                %     constants.r_death = 0.5 * 0.25 * 0.0254;
                %     specified_case = 6;
                %     save_filename = varargin{5};
                
                
                %                 if rem(n_submitted, 450) == 0
                %                     if n_submitted > 0
                %                         pause(30*60);
                %                         % wait 30 minutes to avoid filling qsub que
                %                     end
                %                 end

                [stat, n_l_s] = system('qstat | grep -c ^');
                n_lines = str2double(n_l_s);
                
                if n_lines > 400
                    pause(20*60)
                end
                
                
                error_filename = ['qsubrun' num2str(m) '.error'];
                % check if the error file exists
                if exist(error_filename,'file')
                    
                    % if it does, see if it's size = 0
                    % if it does, then it's already been run with no error
                    s = dir(error_filename);
                    if s.bytes == 0
                        % no error when it was already run
                        run_this_one = 0;
                    else
                        % some error occured. run again
                        run_this_one = 1;
                    end
                else
                    % error file doesn't exist -> run it
                    run_this_one = 1;
                end
                
                
                
                if run_this_one
                    
                    save_filename = ['bubble_sim_data_number_' num2str(m) '.mat'];
                    
                    %             bubble_growth(C_rdot(i), C_nuc_rate(j), C_death_rate(k), C_erf_death_rate(n), save_filename);
                    system(['./code_for_qsub_separate_jobs.sh ' num2str(m) ' ' ...
                        num2str(C_rdot(i)) ' ' num2str(C_nuc_rate(j)) ' ' ...
                        num2str(C_death_rate(k)) ' ' ...
                        num2str(C_erf_death_rate(n)) ' ' save_filename]);
                    %             pause
                    
                    n_submitted = n_submitted + 1;
                end
            end
            
            
        end
    end
end