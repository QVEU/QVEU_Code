# Unix Cheat Sheet
### _Adapted from http://www.mathcs.emory.edu/~valerie/courses/fall10/155/resources/unix_cheatsheet.html (Valerie Summet, Emory)_
## Help on any Unix command.
`man {command}` Retrieve manual entry for command. For example, type `man rm` to read the manual for the `rm` command.

`whatis {command}` Give short description of command.

## _List_ functions
`ls {path}`	It's ok to combine attributes, eg ls -laF gets a long listing of all files with types.

`ls {path_1} {path_2}`	List both {path_1} and {path_2}.

`ls -l {path}` Show all ("long form") information on files. 

`ls -a {path}`	Show all files, including important .dot files that are otherwise hidden.

`ls -F {path}`	Show type of each file. "/" = directory, "*" = executable.

`ls -R {path}`	Recursive listing, with all subdirs.

`ls {path}` | more	Show listing one screen at a time.

## Change to directory
`cd {dirname}`:	__Note:__ There must be a space between.

`cd` or `cd ~`:	Go back to home directory, useful if you're lost.

`cd ..`	Go back up one level in the directory structure, to the parent directory.

## Make a new directory
`mkdir {dirname}`	

## Remove a directory
`rmdir {dirname}`	Only works if {dirname} is empty.

`rm -r {dirname}`	Remove all files and subdirs. __Careful!__

## Print working directory
`pwd`	Show where you are as full path. Useful if you're lost or exploring.

## Copy a file or directory
`cp {file1} {file2}`

`cp -r {dir1} {dir2}`	Recursive, copy directory and all subdirs.

`cat {newfile} >> {oldfile}`	Append newfile to end of oldfile.

## Move (or rename) a file
`mv {oldfile} {newfile}`	Moving a file and renaming it are the same thing.

`mv {oldname} {newname}`	
## Delete a file
`rm {filespec}`	You can use `?` and `*` wildcards to remove files with common names. `?` is any character; `*` is any string of characters.


## View a text file
`more {filename}`	View file one screen at a time.

`less {filename}`	Like more, with extra features.

`cat {filename}`	View file, but it scrolls.

`cat {filename} | more`	View file one screen at a time.

## Edit a text file.
`gedit {filename}`	Basic text editor

`vim {filename}`

`open {filename)`

## Create a text file.
`cat > {filename}`	Enter your text (multiple lines with enter are ok) and press control-d to save.

`gedit {filename}`	Create some text and save it.

`touch {filename}`  Create a file in memory

## Compare two files
`diff {file1} {file2}` Show the differences.

`sdiff {file1} {file2}`	Show files side by side.

## Other text commands
`grep '{pattern}' {file}`	Find regular expression in file.

`spell {file}`	Display misspelled words.

`wc {file}`	Count words in file.

`wc -l {file}`	Count the number of lines in a file.

## Make an Alias
`alias {name}='{command}'` Put the command in 'single quotes'. More useful in your .bashrc file.

## Wildcards and Shortcuts

`*`	Match any string of characters, eg page* gets page1, page10, and page.txt.

`?`	Match any single character, eg page? gets page1 and page2, but not page10.

`[...]`	Match any characters in a range, eg page[1-3] gets page1, page2, and page3.

`~`	Short for your home directory, eg cd ~ will take you home, and rm -r ~ will destroy it.

`.`	The current directory.

`..`	One directory up the tree, eg ls ...

## Pipes and Redirection	(You pipe a command to another command, and redirect it to a file.)
`{command} > {file}`	Redirect output to a file, eg ls > list.txt writes directory to file.

`{command} >> {file}`	Append output to an existing file, eg cat update >> archive adds update to end of archive.

`{command} < {file}`	Get input from a file, eg sort < file.txt

`{command} < {file1} > {file2}`	Get input from file1, and write to file2, eg sort < old.txt > new.txt sorts old.txt and saves as new.txt.

`{command} | {command}`	Pipe one command to another, eg ls | more gets directory and sends it to more to show it one page at a time.

## System info
`date`	Show date and time.

`df`	Check system disk capacity.

`du`	Check your disk usage and show bytes in each directory.

`du -h`	Check your disk usage in a human readable format.

`printenv`	Show all environmental variables.

`uptime`	Find out system load.

`w`	Who's online and what are they doing?

`top`	Real time processor and memory usage

## Unix Directory Format
Long listings (ls -l) have this format:

```
    - file
    d directory,                                            * executable
    ^   symbolic links (?)  file size (bytes)   file name   / directory
    ^           ^               ^                  ^        ^
    drwxr-xr-x 11 valerie    16296  Mar  7 23:25 public_html/
    -rw-r--r--  1 valerie      256  Mar  8 23:42 index.html
                                            ^
     ^^^        user permission  (rwx)      date and time last modified
        ^^^     group permission (rwx)
           ^^^  world permission (rwx)
```
