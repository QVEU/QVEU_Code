
# QVEU Skyline Users Guide

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
