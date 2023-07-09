# EstuaryLarvaSwim
A simple model of estuarine transport and settlement in planktonic larvae

This is a simple model of transport, swimming and settlement of planktonic larvae in an estuary, built on a sophisticated and very useful model of estuarine flows by Parker MacCready (@parkermac) of the School of Oceanography at the University of Washington. This model is presented in detail in two publications:

>MacCready, P. (2004). Toward a unified theory of tidally-averaged estuarine salinity structure. Estuaries, 27, 561-570. doi:10.1007/BF02907644

>MacCready, P. (2007). Estuarine Adjustment. Journal of Physical Oceanography, 37, 2133-2145. doi:10.1175/JPO3082.1
for the time-dependent version.

The generic expectation of flow in estuaries is that river inputs are fresher and lighter than ocean water. River inputs therefore tend to result in surface layers that flow outwards towards the mouth of the estuary. Heavier, salty ocean water tends to form bottom layers that flow inwards towards the head of the estuary. These tendencies are modified by salinity gradients, tidal exchanges and mixing to produce variable regions of transport either retaining larvae inside the estuary or exporting them out into the open ocean. Retention or export of larvae depends strongly on vertical position, and hence on larval swimming behavior.

MacCready's model is an excellent platform in which to explore the consequences for larval dispersal of starting position, time in the plankton and vertical swimming behavior. This is because MacCready's model uses a clever analysis to capture key elements of freshwater inputs from rivers, variations in tidal oscillations and geography, while still requiring minimal computational power. 

The larval swimming model uses codes generously provided by @parkermac to provide the geophysical context for larval release, transport, swimming and settlement in a number of focal estuaries. @parkermac's codes were originally run in MatLab. Here, they have been slightly modified to run in Octave. 

The larval model assumes planktonic larvae of benthic (bottom-living) adults are released just above the bottom of the estuary, within a specified downstream extent that represents adult habitat. After release, larvae undergo a period of planktonic development in the water column. When they have developed for a sufficient time, they become *competent* (capable of settling and metamorphosing into a juvenile or adult form). When competent, larvae settle and metamorphose at their first encounter with adult habitat. If they fail to reach adult habitat during a specified period, or if they get washed out of the estuary, larvae die. The model assesses the horizontal and vertical distributions of larvae that successfully settled, and those that died, at the end of each simulation run.

In addition to the different estuary geographies and levels of river inputs provided by MacCready's model, the larval model focuses on several larval characteristics:

1. Adult habitat:
   - **X_substrate** is the horizontal extent of the estuary in which larvae are released.
   - **Z_substrate** is the vertical extent of depth in which larvae are released (e.g., [0.99 1.] means larvae are released in the bottom 1% of the water column above adult habitat).
   - **settle_only_within_substrate_flag** determines whether larvae can settle anywhere (0) or must settle only within the adult habitat (0).
2. Timing
   - **release_date** is the time (in units of days) that larvae are released, in relation to the start of the simulation and the onset of increased riverine inputs from Spring runoff.
   - **stage1_duration_days** is the number of days in the early development (*pre-competent*) stage of larval development, during which larvae *cannot* yet settle.
   - **stage2_duration_days** is the number of days in the late development (*competent*) stage of larval development, during which larvae *must* settle.
   - **duration_days** is the number of days to run the simulation, after the release date
3. Behavior
   - **v_particle_stage1** is the up- or down-swimming behavior during stage 1, the pre-competent period, with positive velocity in the downwards direction (oceanographers measure vertical position, *Z*, increasing with depth).
   - **v_particle_stage2** is the up- or down-swimming behavior during stage 1, the competent period.
4. Display
   - **colors** determines the color in which larvae are plotted, superimposed on a filled contour plot of horizontal current velocity (e.g. 'r' for red, 'm' for magenta, 'b' for blue, 'c' for cyan).

To faciliate comparing the effects of variations in larval characteristics, the model can simulate several varieties at the same time that vary in some or all of the above characteristics. Variations are set in the model by arrays of parameters. For example, the horizontal distribution of adult habitat for four variants could be

       X_substrate = [-100 -60
                      -90 -50
                      -80 -40
                      -70 -30];  

In this example, variety 1 has 40 km of adult habitat, ranging from -100 km to -60 km from the head of the estuary. Adult habitats for varieties 2, 3 and 4 are successively displaced 10 km towards the mouth of the estuary.

The model plots the magnitude of tidal exchanges, the volume of river inputs, the horizontal and vertical distributions of larvae, and the velocity and larval positions over the estuary. These plots are updated at an interval determined (in units of seconds) by **plot_interval** (default is every two hours).

The model can generate images of each plot, which can be used to make movies of output, by setting **print_plot_flag** to 1 (default is 0, do not make movie images). Images are placed in a directory or folder specified by the variable **directory**, and are given the prefix specifeid by the variable **run_name**. Generating images slows the simulation and takes up considerable disk space, so it's recommended not to make images except when intending to make movies.


