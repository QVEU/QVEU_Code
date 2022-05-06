# Using _locus_, NIAID's HPC cluster
----------------------------------
_Patrick T. Dolan, Ph.D._
_Quantitative Virology and Evolution Unit_

[under construction]

You will need to gain access to locus by going to:

## 1. Log in to _locus_.
- While on campus network or VPN, enter:
``` bash
> ssh username@ai-submit1.niaid.nih.gov
```
- To access the _head node_, you can enter:
``` bash
> ssh username@locus
```
## 2. Work interactively in _locus_ with `qrsh`
- `qrsh` opens an interactive session where you can work on the command line to develop and run analyses.
``` bash
> qrsh -l h_vmem=16G
```
change the memory requirements, `h_vmem`, as needed. 

- this is a good way to develop scripts that you can eventually turn into `qsubmit` job files.

## 3. Transferring files to the cluster with `rsync`
### - To transfer a ___file___:
``` bash
> rsync -avp ~/path/to/local/file username@ai-submit1.niaid.nih.gov:~/path/to/destination/
```

### - To transfer a ___directory___. __Note__: Do not include terminal `/` or you will transfer the _contents_ of the directory into the destination directory.
``` bash
> rsync -avp ~/path/to/local/directory username@ai-submit1.niaid.nih.gov:~/path/to/destination/
```

### - To __download__ data from the cluster back to your machine, just reverse the addresses. 
``` bash
> rsync -avp username@ai-submit1.niaid.nih.gov:~/path/to/data ~/path/to/local/directory/
```


## 4. Submit Jobs in _locus_ with `qsubmit`
___[under construction]___
```
qsubmit jobfile [arguments]
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
----------------------------- /cm/local/modulefiles ------------------------------
cmd             gcc/6.3.0       module-git      openldap
dot             ipmitool/1.8.18 module-info     shared
freeipmi/1.5.5  lua/5.3.4       null

----------------------------- /cm/shared/modulefiles -----------------------------
acml/gcc/64/5.3.1                      intel-cluster-runtime/ia32/2017.0
acml/gcc/fma4/5.3.1                    intel-cluster-runtime/intel64/2017.0
acml/gcc/mp/64/5.3.1                   intel-cluster-runtime/mic/2017.0

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
```

