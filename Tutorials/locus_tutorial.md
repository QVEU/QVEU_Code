Using _locus_, NIAID's HPC cluster
-----------

1. Log in to _locus_.
- While on campus network or VPN, enter:
``` bash
ssh username@ai-submit1.niaid.nih.gov
```
- To access the head  node, you can enter.
``` bash
ssh username@locus
```
2. Work interactively in _locus_ with `sinteractive`
```sinteractive --mem=8g --cpus-per-task=4 --gres=lscratch:30```
