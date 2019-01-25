#!/bin/bash
#
#        FILE: test-convergence.sh
# DESCRIPTION: Runs convergence tests.
#      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
#

echo "=============================================================="
echo " Running convergence tests for xcompact3d"
echo " Author: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>"
echo "=============================================================="

CWD=$(pwd)

MESHES=( 16 )
SSCHEMES=( 1 )
FAILURES="The following tests failed to run:\n"

echo "Running case:"
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
                NSTEP=$(echo "print(int(0.01 / ${dt}))" | python3)

                # Loop over boundary conditions
                for ncx in "${XBCS[@]}"
                do
                    for ncy in "${YBCS[@]}"
                    do
                        # Set mesh size accoring to BCs
                        if [ "${ncx}" = "00" ]; then
                            mshx=${msh}
                        else
                            let mshx=${msh}+1
                        fi
                        if [ "${ncy}" = "00" ]; then
                            mshy=${msh}
                        else
                            let mshy=${msh}+1
                        fi

                        # Setup working directory
                        cd ${CWD}
                        WORK=s${sscheme}/t${tscheme}/b${ncx}${ncy}/n${msh}/dt${dt}
                        echo "- ${WORK}"
                        mkdir -p ${WORK}
                        cp input.i3d ${WORK}
                        cp incompact3d ${WORK}
                        cp probes.prm ${WORK}
                        cp visu.prm ${WORK}
                        cd ${WORK}

                        # Modify input.i3d and run
                        sed -i "s/nx = .*/nx = ${mshx} /g" input.i3d
                        sed -i "s/ny = .*/ny = ${mshx} /g" input.i3d
                        sed -i "s/dt = .*/dt = ${dt} /g" input.i3d
                        sed -i "s/ilast = .*/ilast = ${NSTEP} /g" input.i3d
                        sed -i "s/iorder = .*/iorder = ${sscheme} /g" input.i3d
                        sed -i "s/itimescheme = .*/itimescheme = ${tscheme} /g" input.i3d
                        sed -i "s/nclx1 = .*/nclx1 = ${ncx:0:1} /g" input.i3d 
                        sed -i "s/nclxn = .*/nclxn = ${ncx:1:1} /g" input.i3d 
                        sed -i "s/ncly1 = .*/ncly1 = ${ncy:0:1} /g" input.i3d 
                        sed -i "s/nclyn = .*/nclyn = ${ncy:1:1} /g" input.i3d 

                        mpiexec -np 4 ./incompact3d > OUTPUT.log
                        if [ "$?" != "0" ]; then
                            FAILURES="${FAILURES}${WORK}\n"
                        fi
                    done
                done
            done
        done
    done
done

echo "--------------------------------------------------------------"
echo -e ${FAILURES}
