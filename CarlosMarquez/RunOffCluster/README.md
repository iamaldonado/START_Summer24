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
The file that we will run in the offline cluster, is `RunAnalyses.C`, which is the following. This macro will help us run our analytic train.

The important part that we must remember is the list of which information will be taken. For in the code that is presented is "ListTEST.txt". However, that list only contains 500 events, which were taken to test everything was working well. 

The next section shows how you can change that part for the lists we create.
* [RunAnalyses.C](RunAnalyses.C)

### Code to generated many jobs
This file is the important part so that too much information can be run. hen I will explain it part by part, for the moment is the following and the correct order to be able to run it without any problem.

* [manyjobs.sh](manyjobs.sh)

Once this file is done, We run it using the following command.
```ruby
	source manyjobs.sh
```

### Code to running many files offline cluster

### Code to gather all output files

## Commands to Use Slurm with Offline Cluster

### sbatch

### scancel

### sinfo

### squeue

### Output Files

### Check the Status 
We can use the commands explained in the previous subsection and a mix of them to be able to check the states of what we are running on the cluster offline. 

When we use `squeue` we can see the job-ID, the partition on which it is running, the file name, the user, the status, the time, node and number of the node on which each file is running. 

```ruby
     JOBID  PARTITION     NAME   USER   ST     TIME    NODES NODELIST(REASON)
    1466553 lustre-te runanaun jcmarqzr  R    1:43:36      1 ncx131
    1466552 lustre-te runanaun jcmarqzr  R    1:45:39      1 ncx130
```
We can also combine the "squeue" and "wc" commands to see exactly how many files are running in the offline cluster.

```ruby
    [in]    squeue | wc

    [out]     3      24     235
```

## Possible Errors

### Status of the Cluster


PS
Not running here as it was too much information and filling up. So we decided to run everything in the ~/scratch2/marquez. I will make corrections later to make it easier.
That is the main reason for writting this macro outside.
