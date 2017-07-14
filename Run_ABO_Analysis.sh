source activate minionvironment

#python src/AnalyzeAbo_Main.py  \
# --reference="/home/ben/Github/abo-analysis/input/A1_reference_Exon7.fasta" \
# --alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
# --output="/home/ben/Github/abo-analysis/results/Exon7PolymorphismResults" \
# --analysis-type="ALLELES"
 
#python src/AnalyzeAbo_Main.py  \
# --reference="/home/ben/Github/abo-analysis/input/A1_reference_Exon6.fasta" \
#--alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
# --output="/home/ben/Github/abo-analysis/results/Exon6PolymorphismResults" \
# --analysis-type="ALLELES"

python src/AnalyzeAbo_Main.py  \
 --reference="/home/ben/Github/abo-analysis/input/A1_reference_Exon6.fasta" \
 --alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
 --output="/home/ben/Github/abo-analysis/results/ReadAlignmentResults33332_Ex6" \
 --analysis-type="READS" \
 --reads="/home/ben/Github/abo-analysis/input/33332_prepared/33332_Pass.fastq" \
 
 python src/AnalyzeAbo_Main.py  \
 --reference="/home/ben/Github/abo-analysis/input/A1_reference_Exon7.fasta" \
 --alleles="/home/ben/Github/abo-analysis/input/ABO_Database.fasta" \
 --output="/home/ben/Github/abo-analysis/results/ReadAlignmentResults33332_Ex7" \
 --analysis-type="READS" \
 --reads="/home/ben/Github/abo-analysis/input/33332_prepared/33332_Pass.fastq" \
 
 
source deactivate

