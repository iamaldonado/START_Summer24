# Running in the Offline Cluster

When we have too much information that we want to analyze, a good solution is to divide our data into different files and run them simultaneously, and finally we gather all the output files. To do this we will have to consider several aspects in order to run our analysis train.

## Information and Create Lists
Let's start to our example, we have a lot of information, so I chose to divided on many lists, folders and run to the offline cluster. Where each list have 1111 lines with 500 events.

> Why do we need a lot of lists?

The main reason for making a large list is because we can’t put too much information in one single file `.txt` So, if we try to put all the information. It will send us an error, where it will tell us that the information does not exist. So, we’ll have to cancel everything and divide it into many lists.

### Create many lists
We put the next line in our Linux terminal and with this we will create each list, with the ranges that we want. For our example, we only need to made 20 lists.

First, we use the commmand `ls` by placing the specific address of our files. This help us to find each file we want if we know the address. These are separated by a series of numbers. It's means, the ranges of information that we want.

```ruby
ls /address_of_our_files/*range-number* > name_of_our_new_list.txt
```

Where we put the range between asterisks, since this allows us to take the rest of the files that start with that value. There we will choose from which range we will take the files to form our lists. Then we use `>` to include this entire range in our file `.txt` which will have the name that we decide.

Also, we have to split these lists because it’s a lot of information and could saturate the cluster and make it crash. So we can use the `split` command to divide into lists with the exact number of code lines that we want.

```ruby
split -l Number_of_lines -d name_of_our_list.txt
```
The file is split into short files based on the number of lines which we want using `-l` option as shown. Also, we can change the split files suffix to numeric by using the `-d` option.

## Our file

### Code to generated many jobs

### Code to running many files offline cluster

### Code to gather all output files

## Commands to Use Slurm with Offline Cluster

### sbatch

### scancel

### sinfo

### squeue

## Possible Errors

### Status of the Cluster


PS
Not running here as it was too much information and filling up. So we decided to run everything in the ~/scratch2/marquez. I will make corrections later to make it easier.
That is the main reason for writting this macro outside.
