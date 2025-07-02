# Snakemake_tutorial
```
$ ls -l
total 12
drwxrwxrwx 1 qbaliu qbaliu 4096 Jul  2 17:11 CODE
-rwxrwxrwx 1 qbaliu qbaliu  685 Jul  2 17:10 HowToReadTheNotes.txt
-rwxrwxrwx 1 qbaliu qbaliu 7204 Jul  2 17:10 Notest.txt
-rwxrwxrwx 1 qbaliu qbaliu   20 Jul  2 16:58 README.md
```

After cloning this repository there is only the _CODE/_ directory and some other files.
Next you must **be in the parent directory of this repo /Snakemake_tutorial/** and do the setup
outlined here: https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#requirements

After this you should have the _snakemake-tutorial/_ directory alongside the _CODE/_ directory.
Inside the _snakemake-tutorial/_ directory there should be the _data/_ directory.
Now, go to the _CODE/_ directory, remove the exising _data_ symlink and create a new
one that points to _/Snakemake_tutorial/snakemake-tutorial/data/_ using
_ln -s <path> <link_name>_
