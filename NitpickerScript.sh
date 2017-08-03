source activate minionvironment


cd /home/ben/Github/nit-picker/src
python nit_picker_main.py \
 --reads="/home/ben/Github/abo-analysis/input/33332_raw" \
 --outputdir="/home/ben/Github/abo-analysis/input/33332_prepared" \
 --sampleid="33332" \
 #--minlen="800" \
 #--maxlen="1000" 
 
source deactivate

