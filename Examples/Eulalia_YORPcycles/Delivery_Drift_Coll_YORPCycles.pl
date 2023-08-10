#!/usr/bin/perl 
#use strict;
use Math::Trig;
use constant PI => 3.14159265358979;


#### Simulation parameters
$tstep = 0.01;  # in Myr.
$tmax = 1500;  #Myr
$toutput = 1000; # Print every N timesteps
$tstepcounter = 0;

#collisional cascade
$RadLargeRem = 0.7937;  #Largest remnant see bottke and morby papers
$RadLargeFrag = 0.7937; #fragment is the start of the SFD
$RadSecLarge = $RadLargeFrag; 
$FamSFDslope = -3.5;

# Controlling N in the simulation
$npt = 1000;
$do_coll = 1;
$SmallLimit = 0.5;  # km don't go below this in the loop below
$MaxNew = 1999;  
$MinSize = $SmallLimit;
$MaxSize = 0.7;


#Family/asteroid properties
$pV = 0.045;
#$ainit = 2.4189; #Polana is at 2.4189, 
$ainit = 2.487; # Eulalia is at 2.487au

#Yarko drift properties
$dadtbennu_measured = 19.0e-4; #au/Myr; Drift rate from Spotto + Bolin's TI adjustment
$dadtbennu = $dadtbennu_measured; 
$bennudiam = 0.492;
$bennua    = 1.12;

## YORP cycle variables
$cyorp = 0.6; # 0.6; #from Erigone in Vok2006
$density = 1.1; #from Bennu/Ryugu - needed for YORP
$gammayorp = 0.5; #see page 201 Bottke 2015

##
$do_stoch_yorp = 1;
$tauyorp = 1.0; #Stocahstic YORP timescale in Myr;; DEPRECATED - set by conditional below

#Collisional reorientation
$do_reorient = 1;
$creor = 1.0;  # c_reorient - see Bottke2015


#set global counters
$toadd=0;

# Two outputs. out.out 
open(OUT,">out.out") || die("Could not open file out.out \n");
    print OUT "0 $npt \n";
for ($i=0;$i<$npt;$i++)  {

    # Diam from intiial family SFD - up to Max size.
    $D[$i]=getD();  #km
    $H[$i] =  2.5*(6.244 - log($pV) - 2*log($D[$i]));

    #Select the pole from unit sphere. Ymult always 1 (old code).
    $polgood[$i]=acos(1-2*rand());
    $polsin[$i]=cos($polgood[$i]);
    $Ymult[$i]=1.0; #deprecated here

    
    #set first loop counters
    $col[$i] = 0;
    $form[$i] = 0;
    $gen[$i] = 1;
    $restime[$i] = 0;
    $a[$i] = $ainit;
    $J72[$i] = 0;
    $spin[$i] = getMaxwellSpin(); #start with Maxwell
    $obliq_deg = rad2deg(acos($polsin[$i]));

    #initialize Stochastic YORP timer - see bottke2015 for numbers
    if ($D[$i] > 3.0) {$tauyorp = 1.0;}
    elsif ($D[$i] > 1.9) {$tauyorp = 0.5;}
    elsif ($D[$i] > 0.9) {$tauyorp = 0.25;}
    else {$tauyorp = 0.1;}
    #timer for stoch_yorp activation: between 0-tauyorp
    $stoch_yorp_timer[$i] = rand()*$tauyorp;

    #initiall set updown
    if (rand() < $gammayorp){
	$updown[$i] = 1.0;
    } else {
	$updown[$i] = 1.0;
    }
   

    print OUT "$i $a[$i] $H[$i] $da $adot $Ymult[$i] $D[$i] $restime[$i] $col[$i] $form[$i] $polsin[$i] $gen[$i] $spin[$i] $updown[$i]\n";
    
}

######### BIG LOOP    ########
for ($t=0;$t<$tmax;$t=$t+$tstep) {
    if ($tstepcounter%$toutput == 0) {
	print OUT "$t $npt \n";
	open(TEMP,">temp.temp") || die("Could not open file temp.temp \n");
    }
    
    #From last tstep - how many dudes are getting added. Then reset.
    $npt = $npt + $toadd;
    $toadd = 0;

    #Loop particle by particle.
    for ($i=0;$i<$npt;$i++) {

	# a<0 is the indicator that it has died.
	if ($a[$i] > 0){

	    #get drift in a.
	    $adot = (sqrt($bennua)/sqrt($a[$i]))*(($bennudiam/$D[$i]))*$dadtbennu;
	    $da =  $tstep*$adot*$polsin[$i]*$Ymult[$i];
	    $lasta=$a[$i];
	    $a[$i]=$a[$i]+($da);

	    #Stochastic YORP timer count
	    $stoch_yorp_timer[$i]=$stoch_yorp_timer[$i]+$tstep;
	    
	    #Check on collisional disruption.
	    $destroy = getColl_Ballouz($tstep,$D[$i]);
	    if ($destroy>0 && $do_coll == 1) {

		#Track THIS particle's index through these routines.
		$origindex=$i;
		# This routine determines how many new ones are formed.
		$newnumb=getNumbNew($D[$i]);

		# Loop over this number of new dudes.
		for($j=($npt+$toadd);$j<($npt+$toadd+$newnumb);$j++){

		    getNewD($j,$D[$i],($j-($npt+$toadd)),$origindex);
		    if (rand() < $gammayorp){
			$updown[$j] = 1.0;
		    } else {
			$updown[$j]=1.0;
		    }
		    print "Adding ",$j," ",$D[$i]," ",$D[$j]," ",($j-$npt)," ",$a[$i]," ",$form[$j]," ",$gen[$j]," ",$spin[$j]," ",$updown[$j],"\n";

		}

		$col[$i] = $t;
		$a[$i]=-8;
		$restime[$i] = $t;

		# If we change npt right here... that will be a pain.
		# So we save that number for later.
		$toadd = $toadd + $newnumb;	
	    }
	}
	#Check for orbital drift into resonances - mark their demise.
	if ($a[$i] < 2.15 && $a[$i] > 0) {
	    $a[$i]=-6;
	    $restime[$i]=$t;
	}

	## Now we need to test: On this timestep (da) did we cross
	##   while going the right direction. Then we dont need J72 switch
	##   we can use $lasta as the tester
	if ($a[$i] < 2.2569 && $lasta > 2.2569 && $a[$i] > 0) {
	    #Now, roll the dice. Bottke 2015 says its a 1/3rd chance
	    if (rand() < 0.33333333) {
		$a[$i]=-7;
		$restime[$i]=$t;
	    }
	}
	#if we are moving outward. Note the slightly diff a
	if ($a[$i] > 2.2559 && $lasta < 2.2559 && $a[$i] > 0) {
	    #Now, roll the dice. Bottke 2015 says its a 1/3rd chance
	    if (rand() < 0.33333333) {
		$a[$i]=-77;
		$restime[$i]=$t;
	    }
	}

	if ($a[$i] > 2.499 && $a[$i] > 0) {
	    $a[$i]=-3;
	    $restime[$i]=$t;
	}
	

	if ($a[$i] > 0){

	    #Change rotation rate due to YORP
	    $domdt = domegadt($spin[$i],$D[$i],$a[$i],$density,$polsin[$i]);
	    $omega = (1.0/($spin[$i]*60.0*60.0));
	    $newomega = ($updown[$i]*$tstep*$domdt)+$omega; #Stoch YORP - updown
	    $newspin = 1.0/(($newomega)*60.0*60.0);
	    $spin[$i]=$newspin;

	    # spin, obliq (-1 to 1) go in, and out comes omega de/dt
	    $dedt = dobliqdt($spin[$i],$D[$i],$a[$i],$density,$polsin[$i]);
	    #$omega = (1.0/($spin[$i]*60.0*60.0));
	    $dobdt = 0.5*(180.0/PI)*($dedt); #deg/Myr 0.2 is Bottke hack
	    
	    $newobliq_deg = rad2deg(acos($polsin[$i])) + $tstep*$dobdt;
	    if ($newobliq_deg > 180.0) {
		$newobliq_deg = 179.99;
	    }
	    if ($newobliq_deg < 0) {
		$newobliq_deg = 0.001;
	    } 
	    $polsin[$i] = cos($newobliq_deg/(180.0/PI));

	    ### Check spins  ######
	    ## faster than 2h. stay between 2-4h. 
	    if ($spin[$i] < 2.0) {
		$spin[$i] = rand()*2.0+2.0;
		if (rand() < $gammayorp){
		    $updown[$i] = 1.0;
		} else {
		    $updown[$i] = -1.0;
		}
		#print "FAST SPIN $i - reset to 2\n";
	    }
	    if ($spin[$i] >= 1000.0){
		$spin[$i] = 1000.0;
		#print "SLOW SPIN $i - reset to 1000 - $polsin[$i] $dedt $dobdt $newobliq_deg\n";
		$chck_reorient = 0;
		if ($do_reorient > 0) {
		    $chck_reorient = check_reorient($spin[$i],$D[$i]);
		}
	    
		if($chck_reorient) {

		    $spin[$i]=getMaxwellSpin();
		    #Select the pole from unit sphere. Ymult always 1 (old code).
		    $polgood[$i]=acos(1-2*rand());
		    $polsin[$i]=cos($polgood[$i]);
		    $Ymult[$i]=1.0;
		    if (rand() < $gammayorp){
			$updown[$i] = 1.0;
		    } else {
			$updown[$i] = -1.0;
		    }

		    # print "REOR $i $spin[$i] $polsin[$i] $updown[$i]\n";
		}
		else {#print "slow but no reorient\n";
		}
		
	    }
	} #close the a>0 for re-orienting/domega/dt

	####### STOCHASTIC YORP
	##In the new stochastic YORP formulation, a different torque solu- tion is chosen every time the timescale sYORPðDÞ is exceeded (see Appendix A).
	#Bottke 2015:  0.5 tau=0.1Myr
	#    1km tau = 0.25Myr
	#    2km tau = 0.5Myr
	if ($a[$i] > 0 && $do_stoch_yorp > 0){
	    #simply flips up/down
	    if ($D[$i] > 3.0) {$tauyorp = 1.0;}
	    elsif ($D[$i] > 1.9) {$tauyorp = 0.5;}
	    elsif ($D[$i] > 0.9) {$tauyorp = 0.25;}
	    else {$tauyorp = 0.1;}
	    #has the clock ticked up on stoachastic yorp?
	    if($stoch_yorp_timer[$i] > $tauyorp) {
		if ($i == 1){
		    print"STOCH ",$t," ",$tauyorp," ",$blah," ",$frand," ",$D[$i]," ",$stoch_yorp_timer[$i]," ";
		}

		if (rand() < $gammayorp){
		    $updown[$i] = 1.0;
		} else {
		    $updown[$i] = -1.0;
		}
		if ($i == 1) {
		    print $updown[$i]," \n";
		}
		#reset the timer
		$stoch_yorp_timer[$i]=0;
	    }
		
	}
	
	
	if ($tstepcounter%$toutput == 0) {
	    print OUT "$i $a[$i] $H[$i] $da $adot $Ymult[$i] $D[$i] $restime[$i] $col[$i] $form[$i] $polsin[$i] $gen[$i] $spin[$i] $updown[$i]\n";
	    
	    print TEMP  "$i $a[$i] $H[$i] $da $adot $Ymult[$i] $D[$i] $restime[$i] $col[$i] $form[$i] $polsin[$i] $gen[$i] $spin[$i] $updown[$i]\n";
	}
    }
    #CLEANUP
    $tstepcounter++;
    close TEMP;

}
close OUT;
    print "FINISH ",$npt,"\n";



#### seeing QstarD_Bennu_Equations_MKS in BallouzBoulder folder.
sub getColl_BallouzMKS {
    my $ts = $_[0]; #timestep in year?
    my $di = $_[1]; #dimaeter in km
    my $mu_g=0.33;

    my $U = 5000.0; # m/s MBA impact speed
    my $rho_target = 2000.0; #kg/m**3
    my $rho_imp = 2000.0;
    my $Rw = 80.0; #meters - Otohime.
    
    my $rt = $di*1000.0/2.0; # into MKS

#    my $first = 
    $QstarD = (3.45e2)*$rt**(-0.47) + (2.72e-1)*$rt**(0.99);

}

sub getColl_Ballouz {
    my $ts = $_[0]; #timestep
    my $di = $_[1]; #dimaeter in km
    my $rtarg = $di*100000.0/2.0; #cm and diameter
    my $rtarg_m = $rtarg/100.0;
    my $colverb = 0;
    
    #print "Entering getColl_Ballouz ",$ts," ",$di," ",$rtarg," ",$rtarg_m,"\n";
    
    ###Target Density - need to document in paper
    my $rho_g=1.19; #g/cm3
    my $rho_s=2.2; #g/cm3

    my $mu_g=0.33;

    ### mu_s, we find that it is 0.47 +/- 0.02.
    my $mu_s=0.47;


    #### Q1_func
    my	$G=6.67e-8;
    my	$K1=0.12;
    my	$Kr=1.2;
    my	$mu=$mu_g;
    my  $Rw=8e3; #Otohime radius
    my $ratio_g = 1.0;
    my $delta = 3.0; #gcc
    my $U_mba= 2.5e5; #cm/ss
    ### Q0_func
    my $b=3.0*$mu_g;
    my $a= (-1.0)*$mu_s;

    ### QstarRecon

    

    my $Q1 = 0.5*(($Kr*$K1**(1.0/3.0)**(3.0*(2.0+$mu)/2.0)))*(($G)**(3.0*$mu/2.0))*((4.0*PI/3.0)**($mu-1.0))*($rho_g**(5.0*$mu/2.0))*($delta**(-$mu))*($ratio_g**(3.0*(2.0+$mu)/2.0));
    my $Q0 = -$Q1*$b*($Rw**($b-$a))/($a*$U_mba**(3*$a+$b));

    my $QstarD = ($Q0*($rtarg**$a)*($U_mba**(2.0-(3.0*$mu_s))))+($Q1*($rtarg**$b)*($U_mba**(2.0-(3.0*$mu_g))));

    my $impactor_radius_mba = (($QstarD*($rho_s/$delta)*(2.0/($U_mba)**2))**(1.0/3.0)*$rtarg)/100.0;

    if ($colverb > 1) {
	print "Rtarg Q1 ",$rtarg," ",$Q1,"\n";
	print "Rtarg Q0 ",$rtarg," ",$Q0,"\n";
	print "Rtarg QstarD ",$rtarg," ",$QstarD,"\n";
	print "Rtarg imp_rad_mba ",$rtarg," ",$impactor_radius_mba,"\n";
    }

    
    my $alpha = -2.564813222553626;
    my $C1=13998520270481.287; ##m2.8
    my $Dip1=$impactor_radius_mba*(10.0**0.1);
    
    my $Pi_mba = 2.9e-24; #per m2 per year

    my $N_mba_law = (($C1*$impactor_radius_mba**(-2.8))-($C1*$Dip1**(-2.8)));
    
    my $np_mba = ($N_mba_law*$Pi_mba*(($impactor_radius_mba)+($rtarg/100.0))**2)/2.0;
    
    #lifetime_mba from simplelifetime
    $lifetime = (1.0/$np_mba/1.0e6);
    #####
    
    if ($colverb > 1) {
	print "N_mba_law ",$impactor_radius_mba," ",$N_mba_law,"\n";
	print "Np_mba ",$rtarg," ",$np_mba,"\n";
    }
    if ($colverb > 0){
	$farinella = 16.8*($rtarg_m)**(1.0/2.0); #Myr and meters
	print "Farinella v Ballouz: ",$rtarg_m," ",$farinella," ",$lifetime,"\n";
	#print "Ron Lifetime ",$rtarg," ",$lifetime,"\n";
    }
    
    if ($lifetime < $tstep) { # where we need to do fancy stats
	# you have a VERY short lifetime
	return 1;
    } elsif ($lifetime > $tstep) {
	$prob=($ts/$lifetime); # prob for destruction during this interval
	# THIS SHOULD INCLUDE A NORMAL DIST
	$f=rand();
	if ($prob > $f) {
	    print "remove ",$i," ",$di," ",$lifetime," ",$prob," ",$f," ",$t,"\n"; 
		return 1;
	} else {
	    return 0;
	}
    } else {
	return 0;
    }

}
    
sub getColl {
    my $ts = $_[0];
    my $di = $_[1];

    # lifetime = 16.8 Myr * sqrt(R)
    $lifetime = 16.8*($di*500.0)**0.5;
    #print "getColl ",$di," ",$lifetime," ",$ts/$lifetime,"\n";
    if ($lifetime < $tstep) { # where we need to do fancy stats
	# HACK
	return 1;
    } elsif ($lifetime > $tstep) {
	$prob=($ts/$lifetime); # prob for destruction during this interval
	# THIS SHOULD INCLUDE A NORMAL DIST
	$f=rand();
	if ($prob > $f) {
	    print "remove ",$i," ",$di," ",$lifetime," ",$prob," ",$f," ",$t,"\n"; 
		return 1;
	} else {
	    return 0;
	}
    } else {
	return 0;
    }    
}


##################################################
# Here we want to be able to designate a SFD
#Durda(pg5): 10000 of 2km  1 20km  Super Catastrophic SLOPE = 
# x=0.30103, y=4   x=1.30103 y=0  SLOPE = -4  CUMULATIVE
#Durda(pg5):
### km? YES
#This function gets some variables from the top global settings.
#MinSize
#MaxSize
#It also gets slope from FamSFDslope

sub getD { 
#    $r = rand()*18 + 1;  # sizes from 1-20 km
    
    $small=$MinSize;
    $big=$MaxSize;
    $q=$FamSFDslope;   ### THIS IS THE SLOPE!!!!
    $a=$small**($q);
    $b=$big**($q);
    $f=rand();
    $rd = ((1-$f)*$a + $f*$b)**(1/$q);
#    print "$rd $f $q $a $b \n";
    if ($rd <= 0) {print "Diam = 0"; exit;}
    return $rd;
}


sub getNumbNew {
    my $initD = $_[0];
    my $k=0;
    my $freshD = 0;
    my $count = 0;
    my $massconserve = 0;
    #Question is size of 2nd largest.
    # Answered below. I use the RadSecLarge* the bro$RadLargeRem*$initD;

    #Bootstrap to start with Largest Remnant
    $freshD=$RadLargeRem*$initD;
    $massconserve = ($massconserve**(3.0) + $freshD**(3.0))**(1.0/3.0);
    if ($freshD < $SmallLimit) {
	print "FreshD (RadLargeRem) is below limit - maxed out new guys at ",$count," ",$initD," ",$freshD," ",$massconserve,"\n";
	return $count;
    }

    ## If the first remnant is large enough.. we keep it and start counting.
    $count++;

    for($k=1;$k<$MaxNew;$k++) {

	$freshD = 10**((log(($k+1))/log(10))/($FamSFDslope) + log($RadLargeFrag*$RadLargeRem*$initD)/log(10));
	if ($freshD > $SmallLimit) {
	    $count++;
	    $massconserve = ($massconserve**(3.0) + $freshD**(3.0))**(1.0/3.0);
	} else {
	    print "FreshD (RadLargeFrag) is below limit - maxed out new guys at ",$count," ",$initD," ",$freshD," ",$massconserve,"\n";
	    return $count;
	}

	if ($massconserve > $initD) {
	    print "massconserve limit met getNumbNew ",$initD," ",$freshD," ",$count," ",$massconserve,"\n";
	    return $count;
	}
    }
    print "getNumbNew ",$initD," ",$freshD," ",$count,"\n";

    return $count;
}

    
sub getNewD {
    my $index = $_[0];
    my $initD = $_[1];
    my $cindex = $_[2];
    my $oi = $_[3];
    # WHAT is $i ?? We need to know current npt and then to reset it.

    # Bottke style
    # Largest Remnant first.
    # Then make the largest fragment and teh slope off of it.
    ## CHECK cindex+1
    if ($cindex == 0) {
	$D[$index] = $RadLargeRem*$initD;
    } else {
	$D[$index] = 10**((log(($cindex+1))/log(10))/($FamSFDslope) + log($RadLargeFrag*$RadLargeRem*$initD)/log(10));
    }
	
    $H[$index] =  2.5*(6.244 - log($pV)/log(10) - 2*log($D[$index])/log(10));
    print "getNewD ",$initD," ",$D[$index]," ",$index," ",$cindex,"\n";

    # Random spin poles
    $polgood[$index]=acos(1-2*rand());
    $polsin[$index]=cos($polgood[$index]);
    $Ymult[$index]=1.0;
    
    # Below is for polarized spin.
    #    if (rand () < 0.5) {$polsin[$i]= -1.0;}
    #    else { $polsin[$i]= 1.0;}
    
    
    $col[$index] = 0;
    $form[$index] = $t;
    $restime[$index] = 0;
    $a[$index] = $a[$oi];
    
    #New spin is random too
    $spin[$index] = getMaxwellSpin();

    #set stochastic yorp timer to 0
    $stoch_yorp_timer[$index]=0;
    
    # It is one generation past its parent.
    $gen[$index]=$gen[$oi]+1;
    #print OUT "$i $a[$i] $H[$i] 0 0 $Ymult[$i] $D[$i] $polsin[$i]\n";

}


# Draw from a rayleigh distribution.
# nothing faster than 2hr or slower than 1000h
sub getNewSpin {
    my $temp;
    $rsig = 8.0;
    while (1) {
	$temp = $rsig*sqrt(-2.0*log(rand()));
	if ($temp > 2.0 && $temp < 1000) {
	    return $temp;
	}
    }
}

## Arg. Draw from Maxwellian.
#Johnk’s algorithm
# Verified!!
sub getMaxwellSpin {
    my $temp;
    $msig = 8.0;
    while (1) {
	$r1=rand();
	$r2=rand();
	$r12=$r1*$r1;
	$w=$r1*$r1 + $r2*$r2;
	if ($w < 1) {
	    $temp = -$msig*( ($r12*log(rand())/$w +log(rand()) ) );
	}
	if ($temp > 2.0 && $temp < 1000.0) {
	    return $temp;
	}
    }
}



# From Bottke 2015 Figure 23.
#units are (1/My/s) - (1/sec) (1/Myr). So freq change per million tyear
# TAKE: spin, diameter, semi-major, density, obliq
#returning domega/dt, (1/s)(1/Myr)
sub domegadt {
    my $tspin = $_[0];
    my $tomega = 1.0/($tspin*60.0*60.0); # converts to 1/s; #current omega
    my $di = $_[1]; #dimaeter in km
    my $ta = $_[2]; #current semimajor axis to scale
    my $tdens = $_[3]; #density to scale
    my $tobliq = $_[4];
    

    my $te = rad2deg(acos($tobliq));
    

    my $temp_dwdt = -4.58554e-06+1.55505e-05*sin(0.0361479*($te-46.6356));
    $scale = $cyorp*(2.5/$tdens)*((2.5/$ta)**2.0)*((2.0/$di)**2.0);

    return $scale*$temp_dwdt; #rad to deg
}


# TAKE: spin, diameter, semi-major, density, obliq
#returning radians/Myr
sub dobliqdt {
    my $tspin = $_[0];
    my $tomega = 1.0/($tspin*60.0*60.0); # converts to 1/s; #current omega
    my $di = $_[1]; #dimaeter in km
    my $ta = $_[2]; #current semimajor axis to scale
    my $tdens = $_[3]; #density to scale
    my $tobliq = $_[4]; # probably comes in as 0-1
    my $uudd = $_[5];
    my $te = rad2deg(acos($tobliq)); #switch to deg to match plot

    my $temp_dedt = 1.2e-05*sin(0.0361479*($te-90.0));
    $scale = $cyorp*(2.5/$tdens)*((2.5/$ta)**2.0)*((2.0/$di)**2.0);

    if ($i==1){
	print "DOBLIQ ",$i," ",$t," ",$te," ",$temp_dedt," ",$scale," ",$tspin," ",$adot," ",$a[$i]," ",$updown[$i]," ",$polsin[$i],"\n";
    }
    #returning de/dt.. so omega is removed.
    # this is in 1/Myr.. or radians/Myr -> and now deg
    return (3.14159/180.0)*$scale*$temp_dedt/$tomega;

}

sub check_reorient {
    my $tspin = $_[0];
    my $tomega = 1.0/($tspin*60.0*60.0); # converts to 1/s; #current omega

    my $di = $_[1]; #dimaeter in km

    # If at 1000h, then test tau_reorient = K*P**-5/6 *D**4/3
    #  K is in Vok2006
    # 84.5 kyr  * (omega/omega_0)**5/6 * (D/D_0)**4/3
    ## where D_0 = 2m and omega_0 = 5h
    ### omega_0 = 5.55e-5
    ### d_0 = 0.002 km
    ## 84.5kyr = 0.0845Myr
    ## This si where in Bottke 2015 they go to a Maxwellian with 8h.
    $om_5h = 5.55e-5;
    my $ttreorient = $creor*0.0845*(($tomega/$om_5h)**(5.0/6.0))*(($di/0.002)**(4.0/3.0));
    
    # the probability of happening in THIS timestep: timestep/tau IN MYR
    my $temptemp = $tstep/$ttreorient;
    my $frand = rand();
    #print "chck_reoriented $t $i $tspin $di $tstep $ttreorient $temptemp $frand\n";

    if ($frand < $temptemp) {
	print "reoriented $t $i $tspin $di $tstep $ttreorient \n";
	return 1;
    }
    else {return 0;}
}


############################################################

sub parse {
    my(@parsed);
    @parsed = split(/\s+/,$_[0]);
    chomp(@parsed);
  #remove whitespace at beginning of line, if any
    while ($parsed[0] eq "") {shift(@parsed)}
    return(@parsed);
}
