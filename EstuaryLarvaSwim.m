%
%   This mfile executes a simulation of transport of biotic or abiotic
%   particles in an estuary, as modeled by Parker MacCready's eta2d model.
%
%   This version of the model focuses on the interactions between particle
%   fate and estuarine flows, which vary seasonally as a function of
%   riverine inputs, biweekly as a function of tidal amplitude cycles, and
%   hourly with tidal ebb and flow.
%
%   In the particle transport model, the focus is on how sinking/swimming
%   behaviors, and timing and location of release, affect position in the
%   water column and ultimately distributions in the sediments.
%
%   Particles are transported by advective flows and turbulence (determined
%   by the eta2d model, and also can have intrinsic movement such as sinking
%   or swimming. Particles' intrinsic movements are only vertical. That is,
%   horizontal components of swimming (if there are any) are assumed negligible
%   compared to horizontal advective flows, which are quite substantial. As in
%   the eta2d model, vertical advective flows are neglected, but vertical
%   transport by turbulence can be significant.
%
%   However, particles can have two stages, which are determined by time since
%   release. Movement characteristics can change between these stages.

%   An example of how this can be used is in larval transport, where
%   early developmental stages might swim up and be unable to settle, while
%   stages might swim down and settle on substrate when they encounter it.
%
%   Note that horizontal position units are km, not m (as in the eta2d model).

%   Danny Grunbaum 8/12/2014


%clear all  %  Clear previous results from memory

randn('state',0)  %  Reset the random number generators
rand('state',0)

warning off   %  Get rid of excessive warnings

%  Figure out whether we are running Matlab or FreeMat
FreeMat_flag = isempty(strfind(version,'R20'));
FreeMat_flag = 0 % Set for octave or matlab

tic
%=====================================================



% Specify which of the scenarios below to simulate

% scenario = 0;
scenario = 1;
% scenario = 3;
% scenario = 2;
% scenario = 4;

switch scenario
    case 0   %  Sedimentation: one particle type

        settle_only_within_substrate_flag = 0;  %  This flag determines whether or not particles can settle only
                                                %  within the range of habitat from which they were released.
                                                %
                                                %  settle_only_within_substrate_flag = 1 ==> yes (settle only with original range)
                                                %  settle_only_within_substrate_flag = 0 ==> no (settle on any benthic substrate)
                                                %
                                                %  For example, a larva might need to return to suitable habitat, but a sediment
                                                %  particle might settle anywhere.

        %   Number of particle of different types

        N_part_type = [128];     %   Number of particles to simulate, assorted by type


        release_date = 30;              %   Date on which particles are released

        stage1_duration_days = [4];      %   Number of days particles are in stage 1 (in which they cannot settle)

        stage2_duration_days = [10];  %   Number of days particles are in stage 2 (in which they must settle)

        v_particle_stage1 = [-500e-6];    %

        v_particle_stage2 = [500e-6];    %


        X_substrate = [-100 -80];     %   Horizontal position from which the particles are released
%         X_substrate = [-100 -20];     %   Horizontal position from which the particles are released


        Z_substrate = [.99 1];        %   Range of depths from which particles are released (as fractions of depth)

        colors = ['r'];     %  Colors for each of the types, from b,g,r,c,m,y,k

    case 1   %  Larval dispersal: four larval types

        settle_only_within_substrate_flag = 1;  %  This flag determines whether or not particles can settle only
                                                %  within the range of habitat from which they were released.
                                                %
                                                %  settle_only_within_substrate_flag = 1 ==> yes (settle only with original range)
                                                %  settle_only_within_substrate_flag = 0 ==> no (settle on any benthic substrate)
                                                %
                                                %  For example, a larva might need to return to suitable habitat, but a sediment
                                                %  particle might settle anywhere.

        %   Number of particle of different types

        N_part_type = [64
                       64
                       64
                       64];     %   Number of particles to simulate, assorted by type


        release_date = 30;              %   Date on which particles are released

        stage1_duration_days = [4
                                3
                                2
                                1];      %   Number of days particles are in stage 1 (in which they cannot settle)

        stage2_duration_days = [10
                                8
                                6
                                4];  %   Number of days particles are in stage 2 (in which they must settle)

        v_particle_stage1 = [-500e-6
                             -500e-6
                             -500e-6
                             -500e-6];    %

        v_particle_stage2 = [500e-6
                             500e-6
                             500e-6
                             500e-6];    %


        X_substrate = [-100 -60
                        -100 -60
                        -100 -60
                        -100 -60];     %   Horizontal position from which the particles are released


        Z_substrate = [.99 1.0
                       .99 1.0
                       .99 1.0
                       .99 1.0];        %   Range of depths from which particles are released (as fractions of depth)

        colors = ['r'
                  'm'
                  'b'
                  'c'];     %  Colors for each of the types, from b,g,r,c,m,y,k


    case 2   %  Larval dispersal: four larval types

        settle_only_within_substrate_flag = 1;  %  This flag determines whether or not particles can settle only
                                                %  within the range of habitat from which they were released.
                                                %
                                                %  settle_only_within_substrate_flag = 1 ==> yes (settle only with original range)
                                                %  settle_only_within_substrate_flag = 0 ==> no (settle on any benthic substrate)
                                                %
                                                %  For example, a larva might need to return to suitable habitat, but a sediment
                                                %  particle might settle anywhere.

        %   Number of particle of different types

        N_part_type = [64
                       64
                       64
                       64];     %   Number of particles to simulate, assorted by type


        release_date = 30;              %   Date on which particles are released

        stage1_duration_days = [4
                                3
                                2
                                1];      %   Number of days particles are in stage 1 (in which they cannot settle)

        stage2_duration_days = [10
                                8
                                6
                                4];  %   Number of days particles are in stage 2 (in which they must settle)

        v_particle_stage1 = [-500e-6
                             -100e-6
                             -150e-6
                             -200e-6];    %  swimming up in pre-comp stage

        v_particle_stage2 = [500e-6
                             100e-6
                             150e-6
                             200e-6];    %  swimming down in comp stage


        X_substrate = [-100 -20
                        -80 -20
                        -60 -20
                        -30 -20];     %   Horizontal position from which the particles are released


        Z_substrate = [.99 1
                       .79 .8
                       .59 .7
                       .39 .4];        %   Range of depths from which particles are released (as fractions of depth)

        colors = ['r'
                  'm'
                  'b'
                  'c'];     %  Colors for each of the types, from b,g,r,c,m,y,k




%=====================================================

    case 3   %  Larval dispersal:

        settle_only_within_substrate_flag = 1;  %  This flag determines whether or not particles can settle only
                                                %  within the range of habitat from which they were released.
                                                %
                                                %  settle_only_within_substrate_flag = 1 ==> yes (settle only with original range)
                                                %  settle_only_within_substrate_flag = 0 ==> no (settle on any benthic substrate)
                                                %
                                                %  For example, a larva might need to return to suitable habitat, but a sediment
                                                %  particle might settle anywhere.

        %   Number of particle of different types

        N_part_type = [64
                       64
                       ];     %   Number of particles to simulate, assorted by type


        release_date = 30;              %   Date on which particles are released0

        stage1_duration_days = [8
                                3
                                ];      %   Number of days particles are in stage 1 (in which they cannot settle)

        stage2_duration_days = [8
                                8
                                ];
                                %   Number of days particles are in stage 2 (in which they must settle)

        v_particle_stage1 = [-500e-6
                             -500e-6
                              ];    %  swimming up in pre-comp stage

        v_particle_stage2 = [500e-6
                             500e-6
                             ];    %  swimming down in comp stage


        X_substrate = [-40 -20
                        -40   -20
                        ];     %   Horizontal position from which the particles are released


        Z_substrate = [.79 .8
                       .79 .8
                        ];        %   Range of depths from which particles are released (as fractions of depth)

        colors = ['b'
                  'r'
                   ];     %  Colors for each of the types, from b,g,r,c,m,y,k


    otherwise

        error('Please define this scenario and rerun...')
        end




%=====================================================

%   Estuary characteristics

%load Willapa_low.inp %-mat
load Willapa_medium.inp %-mat
%load Willapa_high.inp
%load test.out
%load Willapa_high.inp -mat

%  Normally, this should be on:
mixing = 1;    %   Turn on turbulent mixing,
% mixing = 0;    %   Turn off turbulent mixing

%=====================================================

%   Simulation parameters

%  Time is in seconds...

start_time = 24*60*60 * release_date;
%end_time =  24*60*60 * max(release_date+stage1_duration_days+stage2_duration_days);
% end_time = 30 * 24 * 3600;
end_time = start_time + 24*60*60 * 10

dt = 15 * 60;  %  This is the time step in seconds
% dt = 60;

small = 1.e-10;
%=====================================================

%   Output/plotting parameters

% plot_interval = 0.25 * 24 * 3600;
%plot_interval = 12 * 3600;
plot_interval = 2 * 3600;
%plot_interval = 1/2 * 3600;

run_name = 'lb2016_duration_';     % This is the prefix for png images, if any
directory = 'images';  % This is the directory in which the images will be saved

%print_plot_flag = 1;  %  This prints out a png file of the particle distribution
                      %  at each plotting interval. The images can be played as
                      %  a movie e.g. with avidemux

                      %  On machines that don't update graphics when requested,
                      %  this forces the update.
print_plot_flag = 0;  %  It takes time, so if it's not needed usually better to
                      %  turn it off.

%  Set velocity contour levels
clevels = [-2:.1:2];
c_axis = [min(clevels) max(clevels)];
contour_label_interval = 2;

%  Generate a new figure to plot results in, and get a "handle"
%  ("swim_fig") so we can bring it to the front when we need it.
disp('got here')
figure
disp('here too')
swim_fig = gcf;
disp('and here')
% whitebg('k')
% whitebg([.5,.5,.5])
if FreeMat_flag==1
    set(swim_fig,'figsize',[800 700])
else
    set(gcf,'position',[50 150 1100 1000]);
%     set(gcf,'position',[50 150 1100 800]);
end

if print_plot_flag == 1
    patch([0 0 1 1],[0 1 1 0],[.5 .5 .5]);
    text(0.25,0.5,'Set figure size and position, then hit <cr>')
    print('junk.png')
    pause
end

t_plot = start_time - plot_interval;





















%=========================================================================
%   INFRASTRUCTURE -- DON'T CHANGE THIS UNLESS YOU KNOW WHAT YOU ARE DOING
%=========================================================================

%  Calculate some parameters, as in Parker's plot_movie_frame code.
fs1 = 10;
% fs1 = 12;
% fs2 = 14;
lw1 = 5;

td = T_vec/86400;
tsd = Tsave_vec/86400;

t = 0; %   Force loading of the first data entry
t_next_time = -1;
% t_index = 0;
t_index = max(find(tsd<release_date))-1;

omega = 2*pi/(3600*12.42);

update_arrays_tidal3

%   Maximum depth is set by the estuary depth at the upstream (left) end of the
%   domain.  Please don't change it unless you have a good reason!

Zmax = -ZZ(1,1);

%=====================================================

%   Initialize particle positions and other variables

N_part = sum(N_part_type)
N_types = length(N_part_type)
N_type_index = [1;cumsum(N_part_type)+1]
% N_type_index = [1;cumsum(N_part_type(1:end-1))+1;N_part]

Xs = zeros(N_part, 1);
Zs = zeros(N_part, 1);

for i = 1:N_types
    Xs(N_type_index(i):N_type_index(i+1)-1) = X_substrate(i,1) + rand(N_part_type(i),1) .* (X_substrate(i,2)-X_substrate(i,1));
end
if FreeMat_flag==1
    Z_bottom = interplin1(X(1,:),-ZZ(1,:),Xs);
    Z_top = interplin1(X(end,:),-ZZ(end,:),Xs);
else
    Z_bottom = interp1(X(1,:),-ZZ(1,:),Xs);
    Z_top = interp1(X(end,:),-ZZ(end,:),Xs);
end

depth = Z_bottom - Z_top;

for i = 1:N_types
    Zs(N_type_index(i):N_type_index(i+1)-1) =  Z_substrate(i,1) * depth(N_type_index(i):N_type_index(i+1)-1) + ...
                    (Z_substrate(i,2)-Z_substrate(i,1)) * depth(N_type_index(i):N_type_index(i+1)-1) .* rand(N_part_type(i),1);
%     Xs(N_type_index(i):N_type_index(i+1)-1) = X_substrate(i,1) + rand(N_part_type(i), 1) * (X_substrate(i,2)-X_substrate(i,1));
end

Us = zeros(N_part, 1);
Ws = zeros(N_part, 1);
Settled = zeros(N_part, 1);

current_day = 0;
current_plot = 0;

%=====================================================
%   Main loop
%=====================================================


for t = start_time:dt:end_time,

   % Update time variables...
    time_in_days = t / (24*60*60);
    if floor(time_in_days)>current_day
        current_day = floor(time_in_days);
        disp(['Day ',num2str(current_day)])
    end


    %   Update geophysical fields, if they are not current
    update_arrays_tidal3

    %   Calculate bottom depth at the X-position of each particle
    if FreeMat_flag==1
        Z_bottom = interplin1(X(1,:),-ZZ(1,:),Xs);
    else
        Z_bottom = interp1(X(1,:),-ZZ(1,:),Xs);
    end
    Z_bottom(find(Xs<X(1,1))) = -ZZ(1,1);
    Z_bottom(find(Xs>X(1,end))) = -ZZ(1,end);

    %   Calculate horizontal and vertical advection, turbulent transport
    %   and swimming movement rates for each particle

    %  Calculate horizontal advection at each particle position; settled particles don't move
    if FreeMat_flag==1
        Us = 1e-3 * bilin_interp(X(1,:),Z(:,1),U,Xs,-Zs./Z_bottom,0.0000001);
    else
        Us = 1e-3 * interp2(X,Z,U,Xs,-Zs./Z_bottom,'*linear',0);
    end
    Us(find(Settled)) = 0;

    %  Calculate eddy diffusivity at each particle position

    if FreeMat_flag==1
        K_eddy_horizontal = mixing * bilin_interp(X(1,:),Z(:,1),KH,Xs,-Zs./Z_bottom,0.0000001);
        K_eddy_vertical = mixing * bilin_interp(X(1,:),Z(:,1),KS,Xs,-Zs./Z_bottom,0.0000001);
    else
        K_eddy_horizontal = mixing * interp2(X,Z,KH,Xs,-Zs./Z_bottom,'*linear',0);
        K_eddy_vertical = mixing * interp2(X,Z,KS,Xs,-Zs./Z_bottom,'*linear',0);
    end

%K_eddy_horizontal = mixing * bilin_interp(X(1,:),Z(:,1),KH,Xs,-Zs./Z_bottom,'*linear',0);
%K_eddy_vertical = mixing * bilin_interp(X(1,:),Z(:,1),KS,Xs,-Zs./Z_bottom,'*linear',0);

    K_eddy_horizontal(find(Settled)) = 0;
    K_eddy_vertical(find(Settled)) = 0;

    %   Move particles through space -- vertical component

    for i = 1:N_types
        if time_in_days < release_date + stage1_duration_days(i)
            v_particle = v_particle_stage1(i);
        else
            v_particle = v_particle_stage2(i);
        end
        Ws(N_type_index(i):N_type_index(i+1)-1) = v_particle;
    end

    Ws(find(Settled)) = 0;
    Us(find(Settled)) = 0;

    %  Make sure "exported" particles can't move
    Us(find(Xs>X(1,end))) = 0;
    Ws(find(Xs>X(1,end))) = 0;

    %  Update particle positions
    Xs = Xs + 1e-3 * randn(size(Xs)) .* sqrt(2*dt.*K_eddy_horizontal) + dt * Us;
    Zs = Zs + randn(size(Xs)) .* sqrt(2*dt.*K_eddy_vertical) + dt * Ws;

    % This is necessary because of a bug in FreeMat's square root function...
    Xs = real(Xs);
    Zs = real(Zs);

    %   Make sure particles can't fly or dig; settle particles that contact substrate

    Zs(find(Zs <= 0)) = small;   %  This insures particles cannot be above the water line

    %  This checks each particle type to determine which particles settle
    for i = 1:N_types

        %  Check if the present time is within the settlement window of particle type i
        if ((time_in_days>release_date+stage1_duration_days(i))&(time_in_days<=release_date+stage1_duration_days(i)+stage2_duration_days(i)))

            if settle_only_within_substrate_flag == 1
                tmp_settle = find( Zs(N_type_index(i):N_type_index(i+1)-1) > Z_bottom(N_type_index(i):N_type_index(i+1)-1) ...
                             & Xs(N_type_index(i):N_type_index(i+1)-1) > min(X_substrate(i,:)) ...
                             & Xs(N_type_index(i):N_type_index(i+1)-1) < max(X_substrate(i,:)));
            else
                tmp_settle = find( Zs(N_type_index(i):N_type_index(i+1)-1) > Z_bottom(N_type_index(i):N_type_index(i+1)-1));
            end

            Settled(tmp_settle+N_type_index(i)-1) = 1;

        end

%         tmp_settle = tmp_settle_all(find(tmp_settle_all));
%
%         if settle_only_within_substrate_flag == 1
%             tmp_settle = find(Zs > Z_bottom & Xs > min(X_substrate) & Xs < max(X_substrate));
%         else
%             tmp_settle = find(Zs > Z_bottom);
%         end
%
%         if time_in_days > release_date + stage1_duration_days
%             Settled(tmp_settle) = 1;
%         end
    end

    tmp = find( Zs>Z_bottom & Settled==0);   % This insures particles cannot be below the top of the substrate
    Zs(tmp) = Z_bottom(tmp) - small;

    %   Make sure particles can't get advected upstream of the left boundary;
    %   this is treated differently than the downstream boundary, where
    %   particles are allowed to "get stuck".

    Xs(find(Xs<X(1,1))) = X(1,1) + small;


    %   Plot the spatial distribution of particles

    if t >= t_plot

        PlotEstuary5d;  %  This mfile contains the plotting commands

        drawnow       % Force the plots to be updated immediately rather
                      % than waiting for the end of the run. On some machines
                      % this seems not to work (a bug) in which case setting
                      % print_plot_flag = 1 will force the update

        if print_plot_flag == 1
           current_plot = current_plot + 1;
           prefix = [directory,filesep,run_name];
           suffix = '.png';
           if current_plot < 10
               fig_filename = [prefix,'000',num2str(current_plot),suffix]
           elseif current_plot < 100
               fig_filename = [prefix,'00',num2str(current_plot),suffix]
           elseif current_plot < 1000
               fig_filename = [prefix,'0',num2str(current_plot),suffix]
           elseif current_plot < 10000
               fig_filename = [prefix,num2str(current_plot),suffix]
           else
               error('Too many plots!')
           end

           if FreeMat_flag==0
               eval(['print -dpng ',fig_filename])
           else
%                print(fig_filename)
               eval(['print ',fig_filename])
           end

        end

    end

end



%   Summary statistics

disp(' ')
disp('Statistics for all particles')
disp(['N = ',num2str(N_part)])

disp(' ')
disp('      Mean       Std')
disp(['X    ',num2str(mean(Xs)),'   ',num2str(std(Xs))])
disp(['Z    ',num2str(mean(Zs)),'   ',num2str(std(Zs))])
disp(' ')

disp(' ')
disp(' ')
disp('Statistics for particles in estuary')

tmp = find(Xs<X(1,end));

disp(['N = ',num2str(length(tmp))])

disp('      Mean       Std')
disp(['X    ',num2str(mean(Xs(tmp))),'   ',num2str(std(Xs(tmp)))])
disp(['Z    ',num2str(mean(Zs(tmp))),'   ',num2str(std(Zs(tmp)))])
disp(' ')




disp(' ')
disp(' ')
disp('Statistics for settled particles')

tmp = find(Xs<X(1,end));

disp(['N = ',num2str(length(tmp))])

disp('      Mean       Std')
disp(['X    ',num2str(mean(Xs(settled_particles))),'   ',num2str(std(Xs(settled_particles)))])
disp(['Z    ',num2str(mean(Zs(settled_particles))),'   ',num2str(std(Zs(settled_particles)))])
disp(' ')

disp(' ')
disp(' ')
disp('Statistics for settled particles, by type')

for i = 1:N_types

    Xs_tmp = Xs(N_type_index(i):N_type_index(i+1)-1);
    Zs_tmp = Zs(N_type_index(i):N_type_index(i+1)-1);
    Xs_in = Xs_tmp(find(Xs_tmp<X(1,end)));
    Zs_in = Zs_tmp(find(Xs_tmp<X(1,end)));

    disp(['Type ',num2str(i),', N = ',num2str(length(Xs_in))])

    disp('      Mean       Std')
    disp(['X    ',num2str(mean(Xs_in)),'   ',num2str(std(Xs_in))])
    disp(['Z    ',num2str(mean(Zs_in)),'   ',num2str(std(Zs_in))])
    disp(' ')
end

toc








