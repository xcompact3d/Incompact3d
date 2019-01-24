#!/bin/bash
#
#        FILE: test-convergence.sh
# DESCRIPTION: Runs convergence tests.
#      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
#

CWD=$(pwd)

MESHES=( 16 32 64 128 )
SSCHEMES=( 1 2 3 4 )

for msh in "${MESHES[@]}"
do
    for sscheme in "${SSCHEMES[@]}"
    do
        if [ "${msh}" == "128" ]; then
            DELTAT=( 4e-3 2e-3 1e-3 5e-4 )
            TSCHEMES=( 1 2 3 4 )
        else
            DELTAT=( 5e-4 )
            TSCHEMES=( 1 )
        fi
        for tscheme in "${TSCHEMES[@]}"
        do
            if [ "${tscheme}" == "1" ]; then
                XBCS=( "00" "11" "22" "12" "21" )
                YBCS=( "00" "11" "22" "12" "21" )
            else
                XBCS=( "00" )
                YBCS=( "00" )
            fi
            for dt in "${DELTAT[@]}"
            do
                # Calculate number of steps to t=2.5
                NSTEP=$(echo "print(int(2.5 / ${dt}))" | python3)

                # Loop over boundary conditions
                for ncx in "${XBCS[@]}"
                do
                    for ncy in "${YBCS[@]}"
                    do
                        # Setup working directory
                        cd ${CWD}
                        WORK=s${sscheme}/t${tscheme}/b${ncx}${ncy}/n${msh}/dt${dt}
                        mkdir -p ${WORK}
                        cp input.i3d ${WORK}
                        cp incompact3d ${WORK}
                        cp probes.prm ${WORK}
                        cp visu.prm ${WORK}
                        cd ${WORK}

                        # Modify input.i3d and run
                        sed -i "s/nx = .*/nx = ${msh} /g" input.i3d
                        sed -i "s/ny = .*/ny = ${msh} /g" input.i3d
                        sed -i "s/dt = .*/dt = ${dt} /g" input.i3d
                        sed -i "s/ilast = .*/ilast = ${NSTEP} /g" input.i3d
                        sed -i "s/iorder = .*/iorder = ${sscheme} /g" input.i3d
                        sed -i "s/itimescheme = .*/itimescheme = ${tscheme} /g" input.i3d
                        sed -i "s/nclx1 = .*/nclx1 = ${ncx:0:1} /g" input.i3d 
                        sed -i "s/nclxn = .*/nclx1 = ${ncx:1:1} /g" input.i3d 
                        sed -i "s/ncly1 = .*/ncly1 = ${ncy:0:1} /g" input.i3d 
                        sed -i "s/nclyn = .*/ncly1 = ${ncy:1:1} /g" input.i3d 

                        mpiexec -np 4 ./incompact3d | tee OUTPUT.log
                    done
                done
            done
        done
    done
done
