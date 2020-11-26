import os

## pre-requisits: python3-pip
## Need install Racon, medaka and config patha

install_pacs = ["clustalx","pbdagcon","build-essential","cmake","ncbi-blast+","clustalw","gcc","make","libbz2-dev","zlib1g-dev","libncurses5-dev","libncursesw5-dev","liblzma-dev"]
pip_pacs = ["numpy","statsmodels==0.9.0","HTSeq","pandas","pybedtools==0.7.10","cython","muscle","biopython==1.74","medaka"]

for pac in install_pacs:
    os.system("sudo apt-get install -y %s"%pac)


for pac in pip_pacs:
    os.system("sudo pip3 install %s"%pac)



