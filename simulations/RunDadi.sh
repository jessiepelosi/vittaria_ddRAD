
#!/bin/bash

echo "
##############################################
##       				    ##
## RUNNING OPTIMIZATION AND DADI INFERENCE  ##
##      	   		            ##
##############################################
"

while getopts "i:r:u:l:o:p:m:" opt; do
	case "$opt" in
		(i) INFILE=${OPTARG};;
		(r) ROUNDS=${OPTARG};;
		(u) MU=${OPTARG};;
		(l) L=${OPTARG};;
		(o) OUTPUT=${OPTARG};;
		(p) RUN=${OPTARG};;
		(m) MODEL=${OPTARG};;
	esac
done

echo "

#############################################
##					   ##
##     RUNNING OPTIMIZATION FOR DADI 	   ##
##      				   ##
#############################################
"

if [ $MODEL = "neutral" ]; then
	python RunDadiOptimization.py -i $INFILE -r $ROUNDS -m neutral | tail -n 20 | grep "Params" | awk '{print $3}' | sed 's/\[//g' | sed 's/\,/ /g' | sed 's/\]//g' > params.txt
	N=$(cat params.txt | awk '{print $1}')

elif [ $MODEL = "two_epoch" ]; then
	python RunDadiOptimization.py -i $INFILE -r $ROUNDS -m two_epoch | tail -n 20 | grep "Params" | awk '{print $3, $4}' | sed 's/\[//g' | sed 's/\,/ /g' | sed 's/\]//g' > params.txt
	NU=$(cat params.txt | awk '{print $1}')
	T=$(cat params.txt | awk '{print $2}')

elif [ $MODEL = "exponential_growth" ]; then
	python RunDadiOptimization.py -i $INFILE -r $ROUNDS -m exponential_growth | tail -n 20 | grep "Params" | awk '{print $3, $4}' | sed 's/\[//g' | sed 's/\,/ /g' | sed 's/\]//g' > params.txt
	NU=$(cat params.txt | awk '{print $1}')
	T=$(cat params.txt | awk '{print $2}')

elif [ $MODEL = "bottlegrowth" ]; then
	python RunDadiOptimization.py -i $INFILE -r $ROUNDS -m bottlegrowth | tail -n 20 | grep "Params" | awk '{print $3, $4, $5}' | sed 's/\[//g' | sed 's/\,/ /g' | sed 's/\]//g' > params.txt
	NUB=$(cat params.txt | awk '{print $1}')
	NUF=$(cat params.txt | awk '{print $2}')
	T=$(cat params.txt | awk '{print $3}')

elif [ $MODEL = "three_epoch" ]; then
	python RunDadiOptimization.py -i $INFILE -r $ROUNDS -m three_epoch | tail -n 20 | grep "Params" | awk '{print $3, $4, $5, $6}' | sed 's/\[//g' | sed 's/\,/ /g' | sed 's/\]//g' > params.txt
	NUB=$(cat params.txt | awk '{print $1}')
	NUF=$(cat params.txt | awk '{print $2}')
	TB=$(cat params.txt | awk '{print $3}')
	TF=$(cat params.txt | awk '{print $4}')

elif [ $MODEL = "three_epoch_inbreeding" ]; then
	python RunDadiOptimization.py -i $INFILE -r $ROUNDS -m three_epoch_inbreeding | tail -n 20 | grep "Params" | awk '{print $3, $4, $5, $6, $7}' | sed 's/\[//g' | sed 's/\,/ /g' | sed 's/\]//g' > params.txt
	NUB=$(cat params.txt | awk '{print $1}')
	NUF=$(cat params.txt | awk '{print $2}')
	TB=$(cat params.txt | awk '{print $3}')
	TF=$(cat params.txt | awk '{print $4}')
	F=$(cat params.txt | awk '{print $5}')

fi

echo "
#############################################
##					   ##
##    	       RUNNING DADI 		   ##
##					   ##
#############################################
"
if [ $MODEL = "neutral" ]; then
	python RunDadi1D.py -i $INFILE -N $N --mutrate $MU --length $L -n $RUN -o $OUTPUT -m neutral

elif [ $MODEL = "two_epoch" ]; then
	python RunDadi1D.py -i $INFILE --nu $NU --T $T --mutrate $MU --length $L -n $RUN -o $OUTPUT -m two_epoch

elif [ $MODEL = "exponential_growth" ]; then
	python RunDadi1D.py -i $INFILE --nu $NU --T $T --mutrate $MU --length $L -n $RUN -o $OUTPUT -m exponential_growth

elif [ $MODEL = "bottlegrowth" ]; then
	python RunDadi1D.py -i $INFILE --nuB $NUB --nuF $NUF --T $T --mutrate $MU --length $L -n $RUN -o $OUTPUT -m bottlegrowth

elif [ $MODEL = "three_epoch" ]; then
	python RunDadi1D.py -i $INFILE --nuB $NUB --nuF $NUF --TB $TB --TF $TF --mutrate $MU --length $L -n $RUN -o $OUTPUT -m three_epoch

elif [ $MODEL = "three_epoch_inbreeding" ]; then
	python RunDadi1D.py -i $INFILE --nuB $NUB --nuF $NUF --TB $TB --TF $TF --F $F --mutrate $MU --length $L -n $RUN -o $OUTPUT -m three_epoch_inbreeding

fi
