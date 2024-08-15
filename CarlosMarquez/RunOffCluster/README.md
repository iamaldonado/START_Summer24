# Running in the Offline Cluster

When we have too much information that we want to analyze, a good solution is to divide our data into different files and run them simultaneously, and finally we gather all the output files. To do this we will have to consider several aspects in order to run our analysis train.

## Information and Create Lists
Let's start to our example, we have a lot of information, so I chose to divided on many lists, folders and run to the offline cluster. Where each list have 1111 lines with 500 events.

> Why do we need a lot of lists?

The main reason for making a large list is because we can’t put too much information in one single file \textit{.txt}. So, if we try to put all the information. It will send us an error, where it will tell us that the information does not exist. So, we’ll have to cancel everything and divide it into many lists.

```ruby
ls /address_of_our_files/*range-number* > name_of_our_new_list.txt
```

PS
Not running here as it was too much information and filling up. So we decided to run everything in the ~/scratch2/marquez. I will make corrections later to make it easier.
That is the main reason for writting this macro outside.
