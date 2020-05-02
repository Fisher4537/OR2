#!/bin/bash

# from tsplib_list
# for each line
# get first world as file name
# curl on the composed url
# grep EDGE_WEIGHT_TYPE | awk '{print $2}'
# if in 

TSP_FILELIST="tsplib_list.txt"
TSP_URL="http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/"
MIN_DIM="300"

while IFS= read -r line;
do
	FILE_NAME=`echo $line | cut -d " " -f1`

	echo "EDGE_WEIGHT_TYPE: `curl -s -G $TSP_URL$FILE_NAME | grep EDGE_WEIGHT_TYPE | awk '{print $3}'`"
	echo "DIMENSION:`curl -s -G $TSP_URL$FILE_NAME | grep DIMENSION | awk '{print $3}'`..."

	EDGE_WEIGHT_TYPE=`curl -s -G $TSP_URL$FILE_NAME | grep EDGE_WEIGHT_TYPE | awk '{print $3}'`
	EDGE_WEIGHT_TYPE_=`curl -s -G $TSP_URL$FILE_NAME | grep EDGE_WEIGHT_TYPE | awk '{print $2}'`

	if [[ "EUC_2D" == $EDGE_WEIGHT_TYPE ]] || [[ "ATT" == $EDGE_WEIGHT_TYPE ]] || [[ "GEO" == $EDGE_WEIGHT_TYPE_ ]]; then
		if [[ "GEO" == $EDGE_WEIGHT_TYPE_ ]]; then
			DIMENSION=`curl -s -G $TSP_URL$FILE_NAME | grep DIMENSION | awk '{print $2}'`
		else
			DIMENSION=`curl -s -G $TSP_URL$FILE_NAME | grep DIMENSION | awk '{print $3}'`
		fi
		
		if [[ $DIMENSION -gt $MIN_DIM ]]; then
			echo "SUCCESSFULLY $FILE_NAME has $EDGE_WEIGHT_TYPE with $DIMENSION nodes..."
			curl -s -G $TSP_URL$FILE_NAME > "data_heavy/$FILE_NAME"
		else
			echo "FAIL DIMENSION CONDITION"
		fi

	else
		echo "$FILE_NAME has $EDGE_WEIGHT_TYPE NOT SUPPORTED..."
	fi
done < "$TSP_FILELIST"


