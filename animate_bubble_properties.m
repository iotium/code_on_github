
t_plot = linspace(0,1,10*round(t(end)*30))*t(end);

movie_save_filename = 'test_movie.avi';

vidobj = VideoWriter(movie_save_filename);

vidobj.FrameRate = 30;

open(vidobj);

% find sauter mean diameter, void fraction, and number density
% versus space, at time indices of t_plot
SMD = NaN*ones(N_nodes,length(t_plot));
alpha = SMD;
ND = SMD;
% dx_nodes = 1/(N_nodes-1);
x_nodes = linspace(0,1,N_nodes)';
r_nodes = linspace(0,1,N_rw);


fig_handle = figure(66);


for i = 1:length(t_plot)
    [~,i_plot] = min(abs(t - t_plot(i)));
    ind_full = 1:(N_full(i_plot));
    SMD(ind_full,i) = 2*mom(ind_full,4,i_plot)./mom(ind_full,3,i_plot);
    alpha(ind_full,i) = V_bubi(ind_full,i_plot);
    ND(ind_full,i) = mom(ind_full,1,i_plot);
    
%     figure(66)
%     subplot(3,3,[1 2])
%     plot(x_nodes, SMD(:,i), 'k')
%     ylabel('SMD')
%     title(['t = ' num2str(t_plot(i))])
% 
%     subplot(3,3,[4 5])
% 
%     plot(x_nodes, alpha(:,i), 'k')
%     ylabel('void frac')
% 
%     subplot(3,3,[7 8])
% 
%     plot(x_nodes, ND(:,i), 'k')
%     ylabel('num den')
%     
%     subplot(3,3,[3, 6, 9])
%     plot(t,P,'b',[t_plot(i) t_plot(i)],[0 P(i_plot)],'k')



    

end

for i = 1:length(t_plot)
            set(fig_handle,'Color',[1 1 1])

    figure(66)
subplot(2,1,1)
plot(x_nodes, ND(:,i), 'k')
xlabel('Normalized Height')
ylabel('Number Density [#/m^3]')
xlim([0 1])
ylim([min(ND(:)), max(ND(:))])

subplot(2,1,2)
plot(x_nodes, 1e3*SMD(:,i), 'k')
xlabel('Normalized Height')
ylabel('Mean Diameter [mm]')
xlim([0 1])
ylim([1e3*min(SMD(:)), 1e3*max(SMD(:))])
    make_text_big('preserve')
    
    if i == 1
        pause
    end
        
    
    
    frame = getframe(fig_handle);
    
    writeVideo(vidobj, frame);
end

close(vidobj)



