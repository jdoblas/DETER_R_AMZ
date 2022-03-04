#!/bin/bash
BASE_DIR="/usr/local/data/DETER_R_AMZ"
TODAY=$(date  '+%Y-%m-%d')
echo "Running DETER_R_AMZ on "$TODAY
echo "Check log on "$BASE_DIR"/logs/detection_"$TODAY".log"
START_DATE=$1
END_DATE=$2
# run app
source /opt/conda/etc/profile.d/conda.sh \
&& conda activate ee \
&& cd $BASE_DIR/src \
&& nohup python -u main.py "$START_DATE" "$END_DATE"  > "$BASE_DIR/logs/detection_$TODAY.log" &
