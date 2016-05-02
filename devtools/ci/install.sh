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
conda config --add channels nividic

conda install --yes rdkit matplotlib pyqt networkx graphviz pygraphviz pil 

conda create -y -n myenv rdkit matplotlib pyqt networkx graphviz pygraphviz pil 

source activate myenv

