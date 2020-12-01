SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Build index files if they don't exists
pisces --config $SCRIPT_DIR/config.json index

# Submit jobs to compute cluster
pisces --config $SCRIPT_DIR/config.json submit -m $SCRIPT_DIR/metadata.csv --overwrite -i gencode_basic gencode_basic_nomask gencode_basic_exonic gencode_basic_exonic_nomask -p 4 --scratch-dir /lustre/scratch/$(whoami)

# Assemble expression matrices
pisces --config $SCRIPT_DIR/config.json summarize-expression -i gencode_basic -m $SCRIPT_DIR/metadata_exp.csv -r Genotype -b Treatment -c Vehicle -d "~Treatment + Genotype" -f $SCRIPT_DIR/contrasts.csv
pisces --config $SCRIPT_DIR/config.json summarize-expression -i gencode_basic_nomask -n expression_matrix_nomask -m $SCRIPT_DIR/metadata_exp.csv -r Genotype -b Treatment -c Vehicle -d "~Treatment + Genotype" -f $SCRIPT_DIR/contrasts.csv
pisces --config $SCRIPT_DIR/config.json summarize-expression -i gencode_basic_exonic -n expression_matrix_exonic -m $SCRIPT_DIR/metadata_exp.csv -r Genotype -b Treatment -c Vehicle -d "~Treatment + Genotype" -f $SCRIPT_DIR/contrasts.csv
pisces --config $SCRIPT_DIR/config.json summarize-expression -i gencode_basic_exonic_nomask -n expression_matrix_exonic_nomask -m $SCRIPT_DIR/metadata_exp.csv -r Genotype -b Treatment -c Vehicle -d "~Treatment + Genotype" -f $SCRIPT_DIR/contrasts.csv

# Make manuscript figures
jupyter nbconvert --to notebook --execute plot_qc_metrics.ipynb

