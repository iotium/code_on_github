% load simulation data and animate plot
movie_save_filename = 'bubble_data_movie.avi';

vidobj = VideoWriter(movie_save_filename);

vidobj.FrameRate = 30;

open(vidobj);

load bubble_sim_data

t_plot = linspace(0,t(end),1000);

for i = 1:length(t_plot)
    
    [~, j] = min(abs(t - t_plot(i)));
    
    N_bubi(N_full(j)+2:end,j) = 1e-16;
    V_bubi(N_full(j)+2:end,j) = 1e-16;
    
    figure(87)
    clf
    x_node = [1:length(V_bubi(:,j))];
    plot(x_node, N_bubi(:,j), 'k-s', x_node, V_bubi(:,j), 'k-o', x_node, (6*V_bubi(:,j)./N_bubi(:,j)/pi).^(1/3),'k-*')
    set(gca,'yscale','log')
    legend('number density','volume density','SMD')
    title(['t = ' num2str(t_plot(i))])

    frame = getframe(gcf);
    
    writeVideo(vidobj, frame);
    
end
    close(vidobj)
close all
