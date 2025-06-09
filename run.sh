# Make the figures directory if it doesn't exist
mkdir -p figures

# If the combined data doesn't exist, download and prepare it
if [ ! -f data/data.csv ]; then
    cd data
    ./download.sh
    Rscript substitutions.R
    python structural_features.py
    Rscript data.R
    cd ..
fi

# Fit the models to the data
Rscript fit.R

# Run all analyses
Rscript analysis/dms_stats.R
Rscript analysis/aa.R
Rscript analysis/aa_physicochemical.R
Rscript analysis/aa_split.R
Rscript analysis/aa_compare.R
Rscript analysis/predictiveness.R
Rscript analysis/surface_accessibility.R
Rscript analysis/demo.R
