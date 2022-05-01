# Using _locus_, NIAID's HPC cluster
----------------------------------
_Patrick T. Dolan, Ph.D._
_Quantitative Virology and Evolution Unit_

[under construction]

You will need to gain access to locus by going to:

## 1. Log in to _locus_.
- While on campus network or VPN, enter:
``` bash
ssh username@ai-submit1.niaid.nih.gov
```
- To access the head  node, you can enter.
``` bash
ssh username@locus
```
## 2. Work interactively in _locus_ with `sinteractive`
- `sinteractive` opens an interactive session where you can work on the command line to develop and run analyses.
```
sinteractive --mem=8g --cpus-per-task=4 --gres=lscratch:30
```

- this is a good way to develop scripts that you can eventually turn into `qsubmit` job files.

## 3. Transferring files to the cluster with `rsync`
- To transfer a file:
```
rsync -avp ~/path/to/local/file username@ai-submit1.niaid.nih.gov:~/path/to/destination/
```
- To transfer a directory. __Note__: Do not include terminal `/` or you will transfer the _contents_ of the directory into the destination directory.
```
rsync -avp ~/path/to/local/directory username@ai-submit1.niaid.nih.gov:~/path/to/destination/
```
## 4. Submit Jobs in _locus_ with `qsubmit`:
TBD.
