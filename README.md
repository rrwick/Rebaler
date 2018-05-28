# Rebaler

Rebaler is a program for conducting reference-based assemblies using long reads. It relies mainly on [minimap2](https://github.com/lh3/minimap2) for alignment and [Racon](https://github.com/isovic/racon) for making consensus sequences.

I made Rebaler for bacterial genomes (specifically for the task of [testing basecallers](https://github.com/rrwick/Basecalling-comparison)). It should in principle work for non-bacterial genomes as well, but I haven't tested it.

* __Reasons to use Rebaler__:
    * You want to reassemble/polish an assembly with long reads, using a reference assembly to guide the large-scale structure.
    * You don't want the reference assembly's basecalls to influence the basecalls in your new assembly.
* __Reasons to _not_ use Rebaler__:
   * You expect there to be structural differences between your reference genome and your sequenced reads. Rebaler is built on the assumption that the genome structure is the same, only bases differ.
   * You want to polish an assembly using the same long reads you used to make the assembly. Rebaler might work for this, but simply using Racon on its own is probably better ([read more here](https://github.com/rrwick/Rebaler/issues/3)).

If you have have raw, signal-level data available (`bax.h5` files for PacBio, `fast5` files for Nanopore), then it's a good idea to run a signal-level polisher ([Arrow](https://github.com/PacificBiosciences/GenomicConsensus) or [Nanopolish](https://github.com/jts/nanopolish)) after Rebaler.





## Requirements

Rebaler runs on [Python 3.4+](https://www.python.org/) and uses [Biopython](http://biopython.org/).

It also assumes that `minimap2` and `racon` executables are available in your `PATH`. If you can open a terminal and run those commands, you're good to go.

Note that as of v0.1.1, Rebaler requires [Racon](https://github.com/isovic/racon) v1.0 or later. Check this with `racon --version`. If that command gives you an error, your version of Racon is too old and you need a new one.



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
    * If you get a strange 'can't combine user with prefix' error, read [this](http://stackoverflow.com/questions/4495120).
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Rebaler`
* Install with pip (from GitHub): `pip3 install git+https://github.com/rrwick/Rebaler.git`


### Run without installation

Rebaler can be run directly from its repository by using the `rebaler-runner.py` script:

```bash
git clone https://github.com/rrwick/Rebaler.git
Rebaler/rebaler-runner.py -h
```



## Usage

Rebaler is simple to use – give it a reference and reads, and it will output its assembly to stdout when done:
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
4) Replace the reference sequence with corresponding read fragments to produce an unpolished assembly (like what [miniasm](https://github.com/lh3/miniasm) would make). If parts of the reference had no read coverage, the original reference sequence will be left in place.
5) Conduct multiple rounds of Racon polishing with all reads to produce the best possible consensus sequence.



## Circular contigs

If the reference is made of circular contigs (as is the norm for bacterial genomes), Rebaler can take this into account during Racon polishing. Specifically, it will 'rotate' the contigs (change the starting position) between polishing rounds to ensure that all parts of the genome are well polished, including the ends.

To indicate that a reference contig is circular, it should have `circular=true` in its fasta header<br>
For example: `>chromosome length=5138942 circular=true`



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
