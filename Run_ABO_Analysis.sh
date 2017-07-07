source activate minionvironment

python src/AnalyzeAbo_Main.py  \
 --reference="/home/ben/Github/abo-analysis/input/A1_reference_Exon7.fasta" \
 --alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
 --output="/home/ben/Github/abo-analysis/results/Exon7PolymorphismResults" \
 --alleles-a="/home/ben/Github/abo-analysis/input/ABO_Database_A.fasta" \
 --alleles-b="/home/ben/Github/abo-analysis/input/ABO_Database_B.fasta" \
 --alleles-o="/home/ben/Github/abo-analysis/input/ABO_Database_O.fasta" \
 
 
source deactivate

