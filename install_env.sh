conda create -n syn_win -y
conda activate syn_win
conda install -c bioconda -c conda-forge python diamond -y
pip install biopython pandas matplotlib seaborn plotly kaleido
