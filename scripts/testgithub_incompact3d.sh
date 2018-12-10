#!/bin/bash
#
### Ask what should be run:-
### - 1 Get the bitbucket repository
### - 2 From an existing bitbucket repository, create all the tests
### - 3 Run on the test which existing folders
#
###
#
echo "*****************************************************************"
echo "      Script to test the various cases of the distribution."
echo "      It starts from the latest 2.0 version of the"
echo "      github.com/xcompact3d/Incompact3d repository, which is"
echo "      downloaded."
echo "*****************************************************************"
#
###
#
rm -rf Incompact3d
rm -rf BENCHMARK_SUITE
rm -rf .output_incompact3d
#
echo
#
mkdir .output_incompact3d
#
read -p "Do you want to use github or local repository? (git/local): " Repo
if [ $Repo == "local" ]
then
		git clone ${I3D_HOME} Incompact3d
else
		git clone https://github.com/xcompact3d/Incompact3d
fi
#
cd Incompact3d
#
echo
read -p "Which branch would you like to use? (default: master): " Branch
if [ -z $Branch ]
then
		echo "Checking out branch " $Branch
		git checkout ${Branch}
fi
#
git status
#
echo
#
git log >& ../OUTPUT_git_log.log
#
git rev-list --all >& ../OUTPUT_git_rev-list.log
#
echo
#
export RevisionNumber=`git rev-list --count HEAD`
#
echo "The latest revision is the $RevisionNumber th one."
#
echo
#
read -p "Do you want to compare with another revision? (Yes/No): " RevisionYesNo
#
if [ $RevisionYesNo == "Yes" ]
    then
    read -p "Revision code for comparison : " OldRevisionCode
#
    OldRevisionLineNumber="$(grep -n "$OldRevisionCode" ../OUTPUT_git_rev-list.log | head -n 1 | cut -d: -f1)"
#
    OldRevisionLineNumber=$(($RevisionNumber - $OldRevisionLineNumber))
    OldRevisionLineNumber=$(($OldRevisionLineNumber + 1))
#
    echo "The revision to compare to is the $OldRevisionLineNumber th."
#
elif [ $RevisionYesNo == "No" ]
    then
    echo "The script will be run for the latest revision number only."
fi
#
echo
#
read -p "Compiler (gcc or intel): " compiler
read -p "Number of processors for the tests: " nprocs
read -p "Number of rows for the tests: " prow
read -p "Number of cols for the tests: " pcol
read -p "Number of time steps: " timesteps
echo
#
### Scan all the BC*.prm files
#
export listBCfiles=`ls examples/*/BC*prm`
#
echo
#
mkdir ../BENCHMARK_SUITE
mkdir ../BENCHMARK_SUITE/$RevisionNumber
#
if [ $RevisionYesNo == "Yes" ]
    then
    mkdir ../BENCHMARK_SUITE/$OldRevisionLineNumber
    listrevision[0]="$OldRevisionLineNumber"
    listrevision[1]="$RevisionNumber"
else
    listrevision[0]="$RevisionNumber"
fi
#
for Revision in $listrevision
do
#
    for prmfile in $listBCfiles
    do
#
####### Identify the tests and get their name
#
        echo "*****************************************************************"
#
        echo
#
        BCfile=${prmfile/.prm/}
        NameTestCase=${BCfile##*/}      # Strip path from case name
				NameTestCase=${NameTestCase:3}  # Strip BC- from case name
#
        echo "   The $NameTestCase test case is prepared using the $compiler compiler"
        echo "   with $nprocs processors and $timesteps time-steps."
#
        echo
#
####### Create a directory for each of the cases (identified from the BC*.prm files)
#
        echo "   A new directory is created, by copying the repository."
        echo "   This directory is called $NameTestCase and located under BENCHMARK_SUITE/$Revision"
        echo "   where $Revision is the revision to be dealt with."
#
        export workdir=`pwd`
#
        mkdir $workdir/../BENCHMARK_SUITE/$Revision/$NameTestCase
        mkdir $workdir/../BENCHMARK_SUITE/$Revision/$NameTestCase/src
        mkdir $workdir/../BENCHMARK_SUITE/$Revision/$NameTestCase/decomp2d
#
####### Copy all the files (Makefile, *.f90, *.prm in each directory
#
        echo
        echo "   The required files (*.f90, *.inc, *.prm and Makefile) are copied under BENCHMARK_SUITE/$Revision/$NameTestCase"
        echo "   and the folders required for the simulation to run are created."
#
        cp decomp2d/*.f90 decomp2d/*.inc ../BENCHMARK_SUITE/$Revision/$NameTestCase/decomp2d/.
				cp src/*.f90 ../BENCHMARK_SUITE/$Revision/$NameTestCase/src/.
				cp $prmfile examples/probes.prm examples/post.prm examples/visu.prm Makefile ../BENCHMARK_SUITE/$Revision/$NameTestCase/.
#
        mkdir ../BENCHMARK_SUITE/$Revision/$NameTestCase/out
#
        cd ../BENCHMARK_SUITE/$Revision/$NameTestCase/
#
####### Change the compiler and flow types for what it needed (ask which one)
####### And add the right option for each case:
#
        echo
#
        sed -i -e "s/.*CMP =.*/CMP = $compiler/" Makefile
        export FLOW_TYPE=$NameTestCase
#
####### Compile and link
#
        echo "   Beginning of compiling/Linking"
#
        make >& compile.log
#
        echo "   End of compiling/linking"
#
        if [ -f "incompact3d" ]
        then
            echo "   The executable file incompact3d was created."
        else
            echo "   The file incompact3d was not created. Check what error(s) occur(s) in the file compile.log."
            echo "   The last 25 lines shows:"
            tail -25 compile.log
        fi
#
        echo
#
####### Change the number of time steps
#
        if [ -f "incompact3d" ]
        then
            sed -i -e "s/.*ilast.*/$timesteps    #ilast/" BC-$NameTestCase.prm
            sed -i -e "s/.*p_row.*/$prow   #p_row/" BC-$NameTestCase.prm
            sed -i -e "s/.*p_col.*/$pcol   #p_col/" BC-$NameTestCase.prm
#
####### Run the simulation in parallel on npcrocs
#
            echo "   Start the simulation for $timesteps time-steps"
#
            mpirun -np $nprocs ./incompact3d >& OUTPUT_${NameTestCase}_${nprocs}_${timesteps}.log
#
        else
            echo "   The simulation will not be run."
        fi
#
####### Complete the run
#
        if [ -f "incompact3d" ]
        then
            echo
            echo
            echo "   The simulation of test case $NameTestCase using the $compiler compiler"
            echo "   with $nprocs processors and $timesteps time-steps is completed."
            echo
            echo
        fi
#
        cd $workdir
#
    done
#
done
#
echo
#
echo "*****************************************************************"
echo "               All the simulations are completed."
echo "               Check the output of each of them"
echo "               in the OUTPUT*.log files"
echo "*****************************************************************"
