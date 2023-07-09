# EstuaryLarvaSwim
A simple model of estuarine transport and settlement in planktonic larvae

This is a simple model of transport, swimming and settlement of planktonic larvae in an estuary, built on a sophisticated and very useful model of estuarine flows by Parker MacCready (@parkermac) of the School of Oceanography at the University of Washington. This model is presented in detail in two publications:

>MacCready, P. (2004). Toward a unified theory of tidally-averaged estuarine salinity structure. Estuaries, 27, 561-570. doi:10.1007/BF02907644

>MacCready, P. (2007). Estuarine Adjustment. Journal of Physical Oceanography, 37, 2133-2145. doi:10.1175/JPO3082.1
for the time-dependent version.

MacCready's model is an excellent platform in which to explore the consequences for larval dispersal of starting position, time in the plankton and vertical swimming behavior. This is because MacCready's model uses a clever analysis to capture key elements of freshwater inputs from rivers, variations in tidal oscillations and geography, while still requiring minimal computational power. 

The larval swimming model uses codes generously provided by @parkermac to provide the geophysical context for larval release, transport, swimming and settlement in a number of focal estuaries. @parkermac's codes were originally run in MatLab. Here, they have been slightly modified to run in Octave. 

The larval model assumes planktonic larvae of benthic (bottom-living) adults are released just above the bottom of the estuary, within a specified downstream extent that represents adult habitat. After release, larvae undergo a period of planktonic development in the water column. When they have developed for a sufficient time, they become *competent* (capable of settling and metamorphosing into a juvenile or adult form). When competent, larvae settle and metamorphose at their first encounter with adult habitat. If they fail to reach adult habitat during a specified period, or if they get washed out of the estuary, larvae die. The model assesses the horizontal and vertical distributions of larvae that successfully settled, and those that died, at the end of each simulation run.

In addition to the different estuary geographies and levels of river inputs provided by MacCready's model, the larval model focuses on several larval characteristics:

1. Adult habitat:
   - **X_substrate** is the horizontal extent of the estuary in which larvae are released
   - **Z_substrate** is the vertical extent of depth in which larvae are released (e.g., [0.99 1.] means larvae are released in the bottom 1% of the water column above adult habitat).
   - **settle_only_within_substrate_flag** determines whether larvae can settle anywhere (0) or must settle only within the adult habitat (0)

