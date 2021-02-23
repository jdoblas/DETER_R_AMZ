#!/bin/bash
TODAY=$(date  '+%Y-%m-%d')
echo "Running SAR_EWS on "$TODAY
echo "Check log on "$HOME"/logs/detection_"$TODAY".log"
START_DATE=$1
END_DATE=$2
# run app
source ~/miniconda3/etc/profile.d/conda.sh \
&& conda activate ee \
&& cd $HOME/CODE/DETER_SAR_EWS/sar_ews \
&& nohup python -u main.py "$START_DATE" "$END_DATE"  > "$HOME/logs/detection_$TODAY.log" &
tail -f "$HOME/logs/detection_$TODAY.log"
