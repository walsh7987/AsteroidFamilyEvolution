# AsteroidFamilyEvolution
Monte Carlo code to model dynamical and collisional evolution of inner Main Belt asteroid families

Author: Kevin Walsh
For publication of Nth generation of Bennu and Ryugu papers

# Software included
This is written in Perl. It is a limited, but forgiving scripting langauage that does not require compilation and should be widely useable and numerous operating systems. It will not be fast, but these are not heavy routines.

## perl scripts

running a perl script from the command line (on a mac or linux based system) is simply:

>perl Delivery_Drift_Coll.pl >& LOG &

where standard output is stashed in LOG, which is useful here. Otherwise out.out are temp.temp are created.


### static obliquities - no YORP cycles
Delivery_Drift_Coll.pl

This version of the code uses static obliquity for each of the evolving simulated asteroids. The simulation parameters are documented in the code, but short explanations are here:

#### simulation parameters
This is the total time to simulate and the timestep in which to do it. Timesteps can be large here since there are not short timescale events (as opposed to YORP cycles where spin rates and obliquities can evolve on ~Myr timescales).

#### collisional cascade
The size of the largest remnant and largest fragment. See Bottke et al. 2002 and Morbidelli et al. 2009 (especially the latter) for discussion of the motivation for these parameters. The former is the size that the disurpted object takes relative to the target. The latter is the size of the next largest object, and it is the anchor of the SFD. The SFD power law slope is also the key parameter.

####family and asteroid properties
Asteroid properities that affect Yarkovsky drift rate. Initial location of the parent body family.

####Yarko drift properties
Yarkovsky drift rates with the known values of Bennu from which drift rate is scaled.


# simulation resolution
The number of objects to start with, the small size limit and the maximum new particles that can be added in any given collisional event. The MaxSize is the maximum size of the initial SFD - the size of the initial parent body's largest remnant. The resolution just affects the statistics of the outcome since each indepdnent particle can be impacted with fragments added.





### YORP cycles
Delivery_Drift_Coll_YORPCycles.pl

## external input files required for run

## examples
Two examples are provided.

###Polana_NoYORP
This script is configured to do a low resolution test of the Polana family for 500Myr with 1000 simulated asteroids. This should take <1min on a modern machine. It is probabilistic by nature so the exact result should vary for each time being run, with the values provided here serving as an example outcome.

Here, out.out is a ~7Mb file. Running Stats.csh on it finds the following:

1494 total final particles
Collision:  862
Inward:  115
Inward nu6:  54
Inward J7:2:  61
Outward:  322
Alive:  195
Alive0:  75
AliveColl:  120
Inward0:  41
InwardColl:  13
Inward J7:2:  38
InwardColl J7:2:  23
Crossover 310 2 



## Quick analysis scripts and plotting files (in Demo)
included is a cshell scrript "Stats.csh" that quickly analyzes the output files.




# System requirements
This is a perl script with minimal resource needs with the exception of storage space for high output or high particle number cases.

Perl can be installed for free and run on nearly any computer platform.

# Installation guide

The scripts are standalone to be run in the Perl command line environment.

# Demo
>perl EXECTUABLE >& OUT &


# Instructions for use - adapting runtime parameters for other needs
The key simulation parameters are described above and in the related publications. The key runtime uses - the output - are described here.

## out.out
This file records the state of each particle in each simulation at each timestep. The timesteps are recorded as:
1. Time
2. Npt

Then the states of the particles are recorded:
1. particle num
2. semi-major axis (negative values are codes for removal)
3. H-mag with assumed albedo
4. drift rate
5. timestep change in a
6. in or out drift
7. Diameter (km)
8. residence time
9. last collision
10. when formed
11. obliquity
12. which generation formed

The removal codes (indicated in column 2) are:
-3 is removal by 3:1
-6 is removal by nu6
-7 is removal by the J72M59
-8 is removal by collisional evolution

## temp.temp
This is the same as above, but only for latest timestep.

## standard output
Lots of useful diagnostics are output to standard output and this is typically saved.


# License is the GNU Open License v3
