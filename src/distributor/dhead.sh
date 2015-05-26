#!/bin/bash
#
# Perform heavy data processing by splitting the task to different servers:
#   Transfer the data and "kernel" function to the servers. 
#   Wait for servers to complete.
#   Collect and merge the results.
#   dhead.sh - "distributed head server shell".

# HEADER INSTRUCTIONS
# ================

#printf "The n input arguments for dhead.sh script are: \n"
#printf "[1] LOGIN : login id to the remote servers (assumed it's the same login for every server)\n"
#printf "[2]- PPATH : working directory (will be created if does not exist) on the remote servers; assumed to be the same for each server \n"
#printf "[3..n-4]- IPADDRS : range of ip-addresses of all the servers, assumed they have the same login/psw account \n"
#printf "[n-3]- REMMAT : name of the matlab function (e.g. 'myfunc') that will be copied and launched on remote servers by dremote.sh \n"
#printf "[n-2]- VARMAT : name of the workspace varialbes file (without numering and .mat); these are the variables to copy and load to matlab memory on the remote servers \n"
#printf "[n-1]- SLEEPTIME : integer that indicates number of seconds to pause when waiting for each remote server to complete their computations \n"
#printf "[n] - FRES : name of the local folder where the results will be copied to from the servers \n"


# ARGUMENTS PARSING
# ================

args=("$@")
printf "\nNumber of arguments passed: %d\n" $#
nargs=$#
if [ $nargs -lt 7 ]; then
	echo "ERROR: Number of passed arguments is smaller than required minimum (7)"
	exit 1
fi
nservs=$((nargs-6))
printf "Number of servers: %d\n" $nservs

LOGIN=${args[0]}
printf "The login parameter: %s\n" $LOGIN
PPATH=${args[1]}
printf "The remote working directory: %s\n" $PPATH
i=2
while true; do
	IPADDRS[$((i-2))]=${args[$i]}
	i=$((i+1))
	nservs=$((nservs-1))
	if (( nservs <= 0 )); then
		break
	fi
done
printf "\nIP addresses extracted:\n"
echo ${IPADDRS[@]}

REMMAT=${args[$nargs-4]} # check file existance
test -e $REMMAT.m
if [ $? -ne 0 ]; then
	printf "ERROR: no such file: %s\n" $REMMAT.m
	exit 1
fi
printf "Matlab script file for remote: %s\n" $REMMAT.m

REMSCRIPT="dserver.sh" # check file existance
test -e $REMSCRIPT
if [ $? -ne 0 ]; then
	printf "ERROR: no such file: %s\n" $REMSCRIPT
	exit 1
fi
printf "Remote bash script: %s\n" $REMSCRIPT

VARMAT=${args[$(($nargs-3))]} # check file existance
j=1
i=0
for IPA in ${IPADDRS[@]}; do
	test -e $VARMAT$j.mat
	if [ $? -ne 0 ]; then
		printf "ERROR: no such file: %s\n" $VARMAT$j.mat
		exit 1
	fi
	FDONE[$i]=0
	IFILES[$i]="$VARMAT$j.mat"
	j=$((j+1))
	i=$((i+1))
done
printf "\nVariables-to-load extracted:\n"
echo ${IFILES[@]}

SLEEPTIME=${args[$(($nargs-2))]}
printf "Pause time is set to %i\n" $SLEEPTIME

FRES=${args[$(($nargs-1))]}
mkdir -p $FRES
printf "The folder to collect result files: %s\n" $FRES
printf "\nFinished reading the input arguments\n"

# CONNECT TO REMOTES, SCP FILES, LAUNCH REMOTE BASH
# ================

eval `ssh-agent`
ssh-add
i=0
printf "\nFile transfer and script launching\n"
for IPA in ${IPADDRS[@]}; do
	ssh $LOGIN@$IPA "mkdir -p $PPATH" # create working directory, if necessary
	ssh $LOGIN@$IPA "rm -f $PPATH/*" # clear the working directory from any previous data
	scp $REMSCRIPT $LOGIN@$IPA:$PPATH 
	scp $REMMAT.m $LOGIN@$IPA:$PPATH 
	scp ${IFILES[$i]} $LOGIN@$IPA:$PPATH 
	
	ssh -n -f $LOGIN@$IPA "sh -c 'cd $PPATH; chmod u+x $REMSCRIPT; nohup ./$REMSCRIPT $REMMAT ${IFILES[$i]} > $VARMAT.out 2> $VARMAT.err < /dev/null &'"
	i=$((i+1))
done

# WAIT LOOP FOR MATLAB FUNCTION ON REMOTE TO TERMINATE
# ================

printf "\nWaiting for Matlab scripts to terminate\n"
TLIMIT=50
count=0
tot=0
while [[ $tot -eq 0 ]]; do
	i=0
	for IPA in ${IPADDRS[@]}; do
		if [ ${FDONE[$i]} -eq 0 ]; then
			printf "Connecting to a server %s and checking for files...\n" $IPA
			ssh $LOGIN@${IPADDRS[$i]} "test -e $PPATH/dserver.dn" # check if *.dn file was generated
			if [ $? -eq 0 ]; then
				FDONE[$i]=1
				printf "Server %d (%s) obtained results\n" $i $IPA
				printf "\nCopying the result files\n"
				nohup scp $LOGIN@$IPA:$PPATH/result_${IFILES[$i]} $FRES & 
			else
				printf " not ready, pause.\n"
				sleep $SLEEPTIME
			fi
		fi
		i=$(($i+1))
	done

	count=$(($count+1))
	if [ $count -eq $TLIMIT ]; then
		printf "Time out - check if requested files exist on remotes\n"
		exit 0
	fi

	tot=1
	i=0
	for IPA in ${IPADDRS[@]} ; do
		if [ ${FDONE[$i]} -eq 0 ]; then
			tot=0
		fi
		i=$(($i+1))
	done
done

kill $SSH_AGENT_PID

printf "\nBash script terminated\n"
