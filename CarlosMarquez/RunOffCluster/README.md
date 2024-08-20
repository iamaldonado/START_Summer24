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

The important part that we must remember is the list of which information will be taken. For in the code that is presented is `ListTEST.txt`. However, that list only contains 500 events, which were taken to test everything was working well. 

The next section shows how you can change that part for the lists we create.
* [RunAnalyses.C](RunAnalyses.C)

### Code to generated many jobs
This file is the important part so that too much information can be run. hen I will explain it part by part, for the moment is the following and the correct order to be able to run it without any problem.
* [manyjobs.sh](manyjobs.sh)

Once this file is done, We run it using the following command.
```ruby
	source manyjobs.sh
```

In the first part we are creating a loop with the `for`, this will start at 0 and will go up to 20. Because with this we will be mapping all the lists that we create with the command `ls`. 

```ruby
	#!/bin/sh

	for ((INDEX = 0; INDEX < 20; INDEX++)) // The number may vary depending on the list we have.
	do
```
After we proceed to create folders inside the folders because we divide our lists and also divide those lists. This is because we had too much information and running it in the offline cluster can become saturated because it can take a lot of time and resources. But by doing that we can decrease the time and resources being used.

To split our lists and also divide those lists. We use the `split` command that we explain in this document. The reason for using it here is because it becomes more efficient, otherwise we would have to divide the lists one by one.

```ruby
	mkdir dir${INDEX}

	cp listareq28_${INDEX}.txt dir${INDEX} 

	cd dir${INDEX}

	split -l 556 -d listareq28_${INDEX}.txt

	mkdir xa0 

	mkdir xa1

	cp x00 xa0/x00
	cp x01 xa1/x01
```

The two new lists that were created are: `xa0` and `xa1`. So we will proceed to create new folders and copy in these folders the new lists. This could be done with a `for` loop. Since that can help make more lists and more folders. However, for the moment it's working, since these lists are sufficient.

However, this is something that will be changed later.

The next part is just copying to folders what you need to run the analysis train without any problem. The part that is important to explain is as follows. With this part of the code, we are telling you that file `RunAnalyses.C` the part called `listTEST.txt` changes it to the name of our new lists. For this case they are `xa0` and `xa1`.

Apart from making this change, you are also copying the modified file to the same folders where the other files were copied so that the analysis train can be run.

```ruby
	sed -e "s/listTEST/x00/" RunAnalyses.C > dir${INDEX}/xa0/RunAnalyses.C
	
	sed -e "s/listTEST/x01/" RunAnalyses.C > dir${INDEX}/xa1/RunAnalyses.C
```

The next part is just to send the runanauno file to the offline cluster with the command `sbatch.

```ruby
	sbatch runanauno
```
This file will be explained in the following subsection.


### Code to running many files offline cluster
The first part of the file shows the general data to run the file in the offline cluster. Where we say the amount of memory and the time needed to be able to do the job without any hitches. As well as in which partition of the offline cluster we want to be executed, in this case it is: nica.

* [runanauno](runanauno)

In the next part, that’s where we’ll load our environment variables. Where we also start the mpdroot and everything you need to run the train without any problems.

And finally we run our analysis train, that is, the file `RunAnalyses.C`.

```ruby
	root -l -b -q "RunAnalyses.C"
```
We use the next commands to run root: `-q` for Exit after processing command line macro files, `-l` for Do not show the ROOT banner, and `-b` for Run in batch mode without graphics.

We are now running many files at once in the offline cluster.

### Code to gather all output files
Once everything we sent to the offline cluster finished running without any problems, we can see that we got a lot of output files. To bring them all together into one and be able to gather the information we want, we use the following file.

Here if we use two loops of `for`, the first one is for the first folders that were created, which as we know go from 0 to 20. The second loop of `for` is from the second division of our lists. They were only divided into 2, so this loop goes from 0 to 2.

The next part, we just tell you where to find our output files, even if a file does not exist just ignore it. In the last part, we tell you how to name the output files with all the information.

Once this file is done. We run it using the Following command.

```ruby
    root mezcla.C
```

## Commands to Use Slurm with Offline Cluster

### sbatch

### scancel

### sinfo

### squeue

### Output Files

### Check the Status 
We can use the commands explained in the previous subsection and a mix of them to be able to check the states of what we are running on the cluster offline. 

When we use `squeue` we can see the job-ID, the partition on which it is running, the file name, the user, the status, the time, node and number of the node on which each file is running.

```ruby9
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
