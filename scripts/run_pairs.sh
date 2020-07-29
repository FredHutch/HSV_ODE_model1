#!/bin/sh
cp hsv_sim.in hsv_sim_fitting.in

if [ "$1" != "" ]
then
    count=$1
else
    count=1
fi
while [ $count -gt 0 ]
do
    echo "$count rounds remaining"
    if [ ! -d "saved_$count" ]
    then
	mkdir "saved_$count"
    fi
    echo "Param_mask 24" >>  hsv_sim_fitting.in
    ./hsv_sim  -b -f hsv_sim_fitting.in -c hsv_sim.crit > c_theta.out
    ../scripts/top_c.pl c_theta.out >> hsv_sim_fitting.in
    ../scripts/top_theta.pl c_theta.out >> hsv_sim_fitting.in
    for file in hsv_sim.out_*
    do
	cp $file "saved_$count/c_theta_$file"
    done

    echo "Param_mask 514" >>  hsv_sim_fitting.in
    ./hsv_sim  -b -f hsv_sim_fitting.in -c hsv_sim.crit > betae_vburst.out
    ../scripts/top_vburstrate.pl betae_vburst.out >> hsv_sim_fitting.in
    ../scripts/top_betae.pl betae_vburst.out >> hsv_sim_fitting.in
    for file in hsv_sim.out_*
    do
	cp $file "saved_$count/betae_vburst_$file"
    done

    echo "Param_mask 5" >>  hsv_sim_fitting.in
    ./hsv_sim  -b -f hsv_sim_fitting.in -c hsv_sim.crit > beta_p.out
    ../scripts/top_beta.pl beta_p.out >> hsv_sim_fitting.in
    ../scripts/top_p.pl beta_p.out >> hsv_sim_fitting.in
    for file in hsv_sim.out_*
    do
	cp $file "saved_$count/beta_p_$file"
    done

    echo "Param_mask 2112" >>  hsv_sim_fitting.in
    ./hsv_sim  -b -f hsv_sim_fitting.in -c hsv_sim.crit > r_eclipse.out
    ../scripts/top_r.pl r_eclipse.out >> hsv_sim_fitting.in
    ../scripts/top_eclipse.pl r_eclipse.out >> hsv_sim_fitting.in
    for file in hsv_sim.out_*
    do
	cp $file "saved_$count/r_eclipse_$file"
    done

    #./printpairs

    count=`expr $count - 1`

done
cp hsv_sim_fitting.in hsv_sim.in.saved
