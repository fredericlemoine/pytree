# pytree

This very preliminary code allows to draw unrooted phylogenetic trees in the "figtree" way.

# Prerequisites

## Biopython

### Under linux or MacOSX

```[bash]
pip install biopython
```

### Under MacOSX

If conda is not installed :

```[bash]
curl -L 'http://repo.continuum.io/miniconda/Miniconda-latest-MacOSX-x86_64.sh' -o Miniconda-latest-MacOSX-x86_64.sh
mkdir -p ~/.anaconda/env/
bash Miniconda-latest-MacOSX-x86_64.sh -b -p ~/.anaconda/env/
export PATH=~/.anaconda/env/bin:$PATH;
conda install anaconda-client
```

You can put the following line in the  `.bash_profile`:

```
export PATH=~/.anaconda/env/bin:$PATH;
```

Then:
```[bash]
conda install biopython
```

## Cairo / pyCairo

### Under linux

```[bash]
sudo apt-get install libcairo2
pip install cairocffi
```

### Under macos

```[bash]
conda install -c https://conda.anaconda.org/vgauthier cairo
conda install -c https://conda.anaconda.org/vgauthier py2cairo
```

# Usage

```
unrooted.py -i <inputfile> -o <outputfile> -t <pdf|png> [-W <width: default 800> -H <height: default 800>]'
```

# Example

```[bash]
cat > t.nw 
((2:1,3:1):1,(4:1,5:1):1,(0:1,1:1):1);

./unrooted.py -i t.nw -o t.png -t png -W 200 -H 200
```

This gives the following tree:
![Example tree](doc/t.png)
