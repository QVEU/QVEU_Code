# QVEU Skyline Users Guide
----------------------------------
_Patrick T. Dolan, Ph.D._
_Quantitative Virology and Evolution Unit_
----------------------------------


This document is intended to describe QVEU-specific 'tips and tricks' for using the NIAID cluster, [Skyline](skyline.niaid.nih.gov). It is not a general UNIX guide (for that, see [UNIX tutorial](https://github.com/QVEU/QVEU_Code/new/main/Tutorials/unix_tutorial.md))

## Make a link to the `lab_share` directory

`ln -s <path to the file/folder to be linked> <the path of the link to be created>`

```
cd
ln -s /hpcdata/lvd_qve/ ./lab_share/
```
`cd` will make sure you are at your home directory `~`, then make a link to the `lvd_qve` home, with a link called `lab_share/`. When you login in the future, you can use this symlink to access the lab directory more easily. it also makes it more easily avilable when using Jupyter notebooks. 

## Set your permissions in the `.bashrc` file
Open the `.bashrc` file with a command like:
```
cd
vim .bashrc
```
and add the line:
```
umask 007
```

## General Use

You will need to request access to skyline by going to: [Skyline](skyline.niaid.nih.gov) on VPN or a campus network. 

## 1. Log in to _skyline_.
- While on campus network or VPN, enter:
``` bash
> ssh username@ai-hpcsubmit1.niaid.nih.gov
```

## 2. Work interactively in _skyline_ with `salloc`
TO DO

## 3. Transferring files to the cluster with `rsync`
### - To transfer a ___file___:
``` bash
> rsync -avp ~/path/to/local/directory username@ai-hpcsubmit1.niaid.nih.gov:~/path/to/destination/dir/
> rsync -avp ~/path/to/local/files/ username@ai-hpcsubmit1.niaid.nih.gov:~/path/to/destination/dir/
> rsync -avp ~/path/to/local/file.ext username@ai-hpcsubmit1.niaid.nih.gov:~/path/to/destination/dir/
```

### - To transfer a ___directory___. __Note__: Do not include terminal `/` or you will transfer the _contents_ of the directory into the destination directory.
``` bash
> rsync -avp ~/path/to/local/directory username@ai-submit1.niaid.nih.gov:~/path/to/destination/
```

### - To __download__ data from the cluster back to your machine, just reverse the addresses. 
``` bash
> rsync -avp username@ai-submit1.niaid.nih.gov:~/path/to/data ~/path/to/local/directory/
```


## 4. Submit Jobs in _locus_ with `sbatch`
```
TO DO
```

## 5. Working with ___modules___
### - To list the currently loaded modules:
``` bash
> module list
```
### - To list available modules:
``` bash
> module avail
```
which returns: 
```
------------------------- /data/apps/modulefiles/misc --------------------------
   ccp4i/8.0         deepemhancer/2023-08-23    phenix/1.21-5207
   ccpem/20221108    eman2/2.99.47
   chimera/1.17.3    fsl/6.0.7.8

----------------------- /data/apps/modulefiles/apptainer -----------------------
   cellpose/2.3.2                immcantation/4.4.0     picard/3.1.1
   colabfold/1.5.3-cuda12.2.2    iphop/1.3.3            vep/110
   gatk/4.5.0.0                  model-angelo/1.0.12    virsorter2/2.2.3

------------------------- /data/apps/modulefiles/conda -------------------------
   acemd/3.7.2-python3.10.12        isoquant/3.3.1-python3.10.12
   ambertools/23.3-python3.10.12    tb-profiler/5.0.1-python3.10.12
   genomad/1.7.3-python3.10.12

------------ /data/apps/modulefiles/spack/linux-rocky9-x86_64/Core -------------
   abseil-cpp/20230125.3-p6svymk
   afni/2024-03-05-rzyorhc
   alsa-lib/1.2.3.2-kpvcosh
   anaconda3/2024.02-1-4lxxvud
...

```
and so on. 

### - To load a module, enter:
``` bash
> module load <modulename>
```
for example, to load minimap2:
``` bash
> module load minimap2
> minimap2
Usage: minimap2 [options] <target.fa>|<target.idx> [query.fa] [...]
Options:
  Indexing:
    -H           use homopolymer-compressed k-mer (preferrable for PacBio)
    -k INT       k-mer size (no larger than 28) [15]
    -w INT       minimizer window size [10]
...
```
