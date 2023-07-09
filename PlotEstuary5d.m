%   This version is modified to omit the vertrical and horizontal profiles
%   in the graphical output.

        figure(swim_fig)
        clf

        t_plot = t_plot + plot_interval;
        edges = [0:2:max(-ZZ(1,:))];
        edges2 = [min(min(X)):20:2];

        %  First, plot all particles
        Xs_in = Xs(find(Xs<X(1,end)));
        Zs_in = Zs(find(Xs<X(1,end)));

        Ns_in = hist(Zs_in,edges);
        Ns_in2 = hist(Xs_in,edges2);

        subplot(3,2,3)

        plot(Ns_in,edges,'k-')
        xlabel('Frequency')
        ylabel('Vertical position')
        title(['Time = ',num2str(t/3600/24),' days'])
        ylim([0 max(-ZZ(1,:))])
        axis ij

        subplot(3,2,4)

        plot(edges2,Ns_in2,'k-')
        xlabel('Horizontal position')
        ylabel('Frequency')
        title(['Time = ',num2str(t/3600/24),' days'])
        axis tight
        axis xy

        %  Now add the corresponding plots by type

        for i = 1:N_types

            Xs_tmp = Xs(N_type_index(i):N_type_index(i+1)-1);
            Zs_tmp = Zs(N_type_index(i):N_type_index(i+1)-1);
            Xs_in = Xs_tmp(find(Xs_tmp<X(1,end)));
            Zs_in = Zs_tmp(find(Xs_tmp<X(1,end)));

            Ns_in = hist(Zs_in,edges);
            Ns_in2 = hist(Xs_in,edges2);

            if length(Ns_in)>0
                subplot(3,2,3)
                hold on
                plot(Ns_in,edges,[colors(i),'-'])
                hold off
            end

            if length(Ns_in)>0
                subplot(3,2,4)
                hold on
                plot(edges2,Ns_in2,[colors(i),'-'])
                hold off
            end

        end

        %  Plot the velocity contours, and superimpose the particle positions, types and states

         subplot(3,1,3)

%          set(gca,'Position',[0.13 0.11 0.775 0.441163])
%          set(gca,'Position',[0.13 0.11 0.775 0.441163])

        h = contourf(X,-ZZ,U,clevels);
        %h = contour(X,-ZZ,U,clevels);
        axis ij
        try
            clim(c_axis)
        catch
            caxis(c_axis)
        end
        clabel(h,clevels(1:contour_label_interval:end))
        xlabel('X (km)')
        ylabel('Z (m)')

        hold on
        plot(X(1,:),-ZZ(1,:),'k-')
        plot(X(1,:),-ZZ(end,:),'k-')

        not_settled_particles = find(Settled==0);
        settled_particles = find(Settled==1);

        for i = 1:N_types

            tmp = not_settled_particles(find( not_settled_particles>=N_type_index(i) & not_settled_particles<=N_type_index(i+1)-1));

            if (time_in_days<=release_date+stage1_duration_days(i))
                %  Particle is in stage 1
                plot(Xs(tmp),Zs(tmp),'k^','markersize',5,'markerfacecolor',colors(i))
            elseif ((time_in_days>release_date+stage1_duration_days(i))&(time_in_days<=release_date+stage1_duration_days(i)+stage2_duration_days(i)))
                % Particle is in stage 2, but not settled
                plot(Xs(tmp),Zs(tmp),'kv','markersize',5,'markerfacecolor',colors(i))
            else
                % Particle is beyond stage 2, and can no longer settle
                plot(Xs(tmp),Zs(tmp),'ko','markersize',5,'markerfacecolor',colors(i))
            end

            tmp = settled_particles(find( settled_particles>=N_type_index(i) & settled_particles<=N_type_index(i+1)-1));
            plot(Xs(tmp),Zs(tmp),'ks','markersize',5,'markerfacecolor',colors(i))

        end

        xlim([min(X(1,:)) max(X(1,:))+5])
        hold off

        Xs_in = Xs(find(Xs<X(1,end)));
%         title(['Time = ',num2str(t/3600/24),' days; ',num2str(N_part-length(Xs_in)),'/',num2str(N_part),' particles exported;  ',num2str(sum(Settled)),'/',num2str(N_part),' particles settled'])
        title([num2str(N_part-length(Xs_in)),'/',num2str(N_part),' particles exported;  ',num2str(sum(Settled)),'/',num2str(N_part),' particles settled'])





    % FORCING
    subplot(3,1,1)
    set(gca,'Position',[0.13 0.683837 0.775 0.241163])
    plot(td,Qr_vec/1000,'-b',td,Ut0_vec,'-r',...
        tsd(ii),this_qr/1000,'ob',tsd(ii),this_ut,'or', ...
        'linewidth',lw1);
    set(gca,'fontsize',fs1);
    axis([0 max(tsd) 0 2.1]);
    xlabel('Time (days)','fontsize',fs1);
    title(['Current Time = ',num2str(t/3600/24),' days'])
    %
    [xt,yt] = pmlab(axis,'ur');
    text(xt,yt,'River Flow / 1000 m^{3} s^{-1}','fontsize',fs1,'color','b', ...
        'horizontalalignment','right','fontweight','bold')
    [xt,yt] = pmlab(axis,'ul');
    text(xt,yt,'Tidal Velocity (m s^{-1})','fontsize',fs1,'color','r', ...
        'horizontalalignment','left','fontweight','bold')
