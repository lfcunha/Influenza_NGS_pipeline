# README #

Repository for de novo influenza genome assembly pipeline.
<p>Original pipeline by Harm Van Bakel</p>
<p>Contribution with the python modules by Luis Cunha</p>

## Requirements

The assembly pipeline depends on a number of packages that are availble on minerva as loadable modules. 
In addition it uses scripts from a few other git repositories that must be available in your path. The 'run-flu-assembly' script will take care of loading all the required modules, but you'll have to set up the additionally required git repositories. To install the influenza assembly pipeline on minerva, create or go to the directory where you keep your git repositories and execute the following:

```
#!bash
cd $GIT_REPODIR
$ git clone git@bitbucket.org:hvbakel/influenza-ngs-assembly.git
$ git clone git@bitbucket.org:hvbakel/igb-tools.git
$ git clone git@bitbucket.org:hvbakel/ngs-tools.git
$ git clone git@bitbucket.org:hvbakel/utility.git
```
Then set the following environment variables in your .bashrc:


```
#!bash
export GIT_REPODIR=/path/to/repodir
export PATH="$GIT_REPODIR/influenza-ngs-assembly/bin:$PATH"
```


## Usage

Once the pipeline is installed, you can run assemblies with the following options:

Assemble a new genome:

```
#!bash

run-flu-assembly in=<sample>.fastq.gz assemble
```


Prepare prepare assembly report for a curated genome:

```
#!bash

run-flu-assembly in=<sample>.fastq.gz report
```


Discover low-frequency snv and indel variants with the lofreq package:

```
#!bash

run-flu-assembly in=<sample>.fastq.gz lofreq
```


Annotate a new genome using the NCBI flu annotation tool:

```
#!bash

run-flu-assembly in=<sample>.fastq.gz annotate
```

   
Upload a new flu genome into cripdb:

```
#!bash

run-flu-assembly in=<sample>.fastq.gz upload
```
