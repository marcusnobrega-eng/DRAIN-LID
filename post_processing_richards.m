% Richards Model -Post-Processing
% Developer: Marcus Nobrega
% Goal - Solve 1-D Richards with Known Boundary Conditions
% Next Steps: Implement a hydrologic model and include boundary conditions
% to allow ponding depth
% Moreover, implement pollutant transport model

%% Plotting Results
number_of_plots = 12;
time_result = tfinal*60/number_of_plots; % min
time_plot = 0:time_result*60/dt:(tfinal*60*60/dt); % Time-step index
time_plot(1,1) = 1;
time_plot(1,end) = time_plot(1,end)-1;
labels = strcat(num2str(round(time_plot*dt/3600,2)),'h');
figure(2)
ax1 = subplot(3,1,1);
n_plots = length(time_plot);
colors = linspecer(n_plots,1);
set(gcf,'units','inches','position',[1,1,6.5,7])
for i = 1:n_plots
    txt = ['$t =~$',num2str(round(time_plot(i)*dt/3600,2)),'$~h$'];
    plot(theta(:,time_plot(i)),(-1)*soil_depth,'color',colors(i,:),'LineStyle','-','LineWidth',2,'DisplayName',txt)
    xlabel('$\theta$ [$\mathrm{cm^3.cm^{-3}}$]','interpreter','latex','FontSize',12);
    ylabel('z [cm] - From Top','interpreter','latex','FontSize',12);
    grid on
    hold on
end
l = legend;
l.Interpreter = 'latex';
l.Location = 'northoutside';
l.NumColumns = 5;
hold off

ax2 = subplot(3,1,2);
for i = 1:length(time_plot)
    txt = ['$t =~$',num2str(round(time_plot(i)*dt/3600,2)),'$~h$'];
    plot(head(:,time_plot(i)),(-1)*soil_depth,'color',colors(i,:),'LineStyle','-','LineWidth',2,'DisplayName',txt)
    xlabel('Head [cm]','interpreter','latex','FontSize',12);
    ylabel('z [cm] - From Top','interpreter','latex','FontSize',12);
    grid on
    hold on
end
l = legend;
l.Interpreter = 'latex';
l.Location = 'northoutside';
l.NumColumns = 5;
hold off

ax2 = subplot(3,1,3);
for i = 1:length(time_plot)
    txt = ['$t =~$',num2str(round(time_plot(i)*dt/3600,2)),'$~h$'];
    plot(q(:,time_plot(i))*86400,(-1)*soil_depth,'color',colors(i,:),'LineStyle','-','LineWidth',2,'DisplayName',txt)
    xlabel('Flux [cm/day]','interpreter','latex','FontSize',12);
    ylabel('z [cm] - From Top','interpreter','latex','FontSize',12);
    grid on
    hold on
end
l = legend;
l.Interpreter = 'latex';
l.Location = 'northoutside';
l.NumColumns = 5;
hold off

exportgraphics(gcf,'Results.pdf','ContentType','vector')

%% Animation Theta
% Video
figure(3)
obj = VideoWriter('THETA.avi','Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
time_frame = 5; % min
number_of_frames = tfinal*60/time_frame;
for n=1:1:number_of_frames
    if n == 1
        t = 1;
        pos = 1;
    else
        t=(n-1)*time_frame*60; % min
        pos = t/dt;
    end
    x = linspace(max(theta(:,pos)),min(theta(:,pos)),nz);
    y = 0*soil_depth;
    plot(x,y,'LineWidth',4,'LineStyle','-','Color','k')
    hold on
    plot(theta(:,pos),(-1)*soil_depth,'k','LineWidth',2,'LineStyle','-','Color','blue')
    xlabel('$\theta$ [$\mathrm{cm^3.cm^{-3}}]$','Interpreter','latex');
    ylabel('$z$ [m] - From Top','Interpreter','latex');
    %     ylim([0.98*min(min(wse - y)) max(max(1.01*wse))])
    grid on
    title(['t = ',num2str(round(round(t/3600,2),2)),' [h]'])
    f = getframe(gcf);
    writeVideo(obj,f);
    hold off
end
obj.close();

%% Animation Head
% Video
figure(4)
obj = VideoWriter('HEAD.avi','Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
time_frame = 5; % min
number_of_frames = tfinal*60/time_frame;
for n=1:1:number_of_frames
    if n == 1
        t = 1;
        pos = 1;
    else
        t=(n-1)*time_frame*60; % min
        pos = t/dt;
    end
    x = linspace(max(theta(:,pos)),min(theta(:,pos)),nz);
    y = 0*soil_depth;
    plot(x,y,'LineWidth',4,'LineStyle','-','Color','k')
    hold on
    plot(head(:,pos),(-1)*soil_depth,'k','LineWidth',2,'LineStyle','-','Color','red')
    xlabel('$\mathrm{head}$ [$\mathrm{cm}]$','Interpreter','latex');
    ylabel('$z$ [m] - From Top','Interpreter','latex');
    %     ylim([0.98*min(min(wse - y)) max(max(1.01*wse))])
    grid on
    title(['t = ',num2str(round(round(t/3600,2),2)),' [h]'])
    f = getframe(gcf);
    writeVideo(obj,f);
    hold off
end
obj.close();

%% Animation q
% Video
figure(4)
obj = VideoWriter('flow.avi','Motion JPEG AVI');
obj.Quality = 100;
obj.FrameRate = 20;
open(obj)
time_frame = 5; % min
number_of_frames = tfinal*60/time_frame;
for n=1:1:number_of_frames
    if n == 1
        t = 1;
        pos = 1;
    else
        t=(n-1)*time_frame*60; % min
        pos = t/dt;
    end
    plot(q(:,pos)*86400,(-1)*soil_depth,'k','LineWidth',2,'LineStyle','-','Color','green')
    xlabel('$\mathrm{q}$ [$\mathrm{cm/day}]$','Interpreter','latex');
    ylabel('$z$ [m] - From Top','Interpreter','latex');
    %     ylim([0.98*min(min(wse - y)) max(max(1.01*wse))])
    grid on
    title(['t = ',num2str(round(round(t/3600,2),2)),' [h]'])
    f = getframe(gcf);
    writeVideo(obj,f);
    hold off
end
obj.close();