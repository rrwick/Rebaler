# Rebaler

Rebaler is a program for conducting reference-based assemblies using long reads. It relies mainly on [minimap2](https://github.com/lh3/minimap2) for alignment and [Racon](https://github.com/isovic/racon) for making consensus sequences.

I made Rebaler for bacterial genomes (specifically for the task of [testing basecallers](https://github.com/rrwick/Basecalling-comparison)). It should in principle work for non-bacterial genomes as well, but I haven't tested it.



## Requirements

Rebaler runs on [Python](https://www.python.org/) 3.4+ and uses [Biopython](http://biopython.org/).

It also assumes that `minimap2` and `racon` executables are available in your `PATH`. If you can open a termanal and run those commands, you're good to go.



## Installation

### Install from source

Running the `setup.py` script will install a `rebaler` executable:

```bash
git clone https://github.com/rrwick/Rebaler.git
cd Rebaler
python3 setup.py install
rebaler -h
```

* If the `python3 setup.py install` command complains about permissions, you may need to run it with `sudo`.
* Install just for your user: `python3 setup.py install --user`
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Rebaler`
* Install with pip (from GitHub): `pip3 install git+https://github.com/rrwick/Rebaler.git`


### Build and run without installation

The program can then be executed by directly from the repository directory using the `rebaler-runner.py` script:

```bash
git clone https://github.com/rrwick/Rebaler.git
cd Rebaler
./rebaler-runner.py -h
```



## Usage

Rebaler is simple to use, just give it reads and a reference and it will output its assembly to stdout when done:
```
rebaler reference.fasta reads.fastq.gz > assembly.fasta
```

Progress information will be outputted to stderr.

### Full usage
```
usage: rebaler [-h] [-t THREADS] [--keep] reference reads

Rebaler: reference-based long read assemblies of bacterial genomes

positional arguments:
  reference                      FASTA file of reference assembly
  reads                          FASTA/FASTQ file of long reads

optional arguments:
  -h, --help                     show this help message and exit
  -t THREADS, --threads THREADS  Number of threads to use for alignment and polishing (default: 8)
  --keep                         Do not delete temp directory of intermediate files (default: delete
                                 temp directory)
```



## Method

<p align="center"><img src="misc/diagram.png" alt="method diagram" width="100%"></p>

1) Load in the reference contigs.
2) Use minimap2 to align long reads to the reference.
3) Remove lower quality alignments (judged by length, identity and size of indels) until the reference is just covered. Any given position in the reference should now have a coverage of 1 or 2 (or 0 if the reads failed to cover a spot).
4) Replace the reference sequence by read sequences to produce an unpolished assembly, similar to what [miniasm](https://github.com/lh3/miniasm) would make.
5) Finally, conduct multiple rounds of Racon polishing with all reads to produce the best possible consensus sequence.



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
