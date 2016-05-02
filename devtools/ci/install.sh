MINICONDA=Miniconda2-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b -p $HOME/miniconda
PIP_ARGS="-U"

export PATH=$HOME/miniconda/bin:$PATH

sudo apt-get update

conda update --yes conda
conda config --add channels http://conda.binstar.org/nividic

conda install --yes rdkit matplotlib networkx graphviz pygraphviz pil pyqt

conda create -y -n myenv python=$PYTHON_VERSION rdkit matplotlib networkx graphviz pygraphviz pil pyqt

source activate myenv


