## Overview

This is a generator for chess endgame databases ("tablebases")
for up to 8 pieces on reasonable hardware.

At the moment, only pawnless tablebases and tablebases with exactly one pawn
are supported.

The generated tablebase files have the same file extensions as regular Syzygy
tablebase files but are not compatible with current chess software (with the
exception of [https://github.com/syzygy1/probetool/tree/new_format](https://github.com/syzygy1/probetool/tree/new_format)). Existing
chess software will likely not crash but will instead output an error or warning
when trying to access a file in the new format (since the "magic number" in the
file header is different), or simply ignore the file.

> **Warning:** Generating an 8-piece tablebase file with this generator will
> probably take many weeks, if not months. The compressed tablebase files will
> be huge and essentially impossible to distribute over the internet (the total
> size of the complete set of 8-piece WDL and DTZ tablebase files will likely
> be between 1500 and 2000 TB). I do not really expect people to try this until
> data storage technology significantly improves, but we will see.

> **Warning:** This generator is still a work in progress. I do not guarantee
> that the current format is final or even that the generated tables are
> correct. In particular, updates may temporarily break the generator. That
> said, if a change in the format turns out to be necessary after a number of
> 8-piece tables have already been generated, I will look into keeping the
> probing code backward compatible and/or writing a tool to convert generated
> tables to the updated format.

## Hardware Requirements

Only Linux systems with an x86-64 CPU supporting efficient PDEP and PEXT
instructions are supported (AMD Zen 3+ and Intel Haswell or newer).

For a pawnless 8-piece tablebase with no identical pieces of the same
color, such as `KQRBvKQRN`, the generator requires a PC with 128 GB of RAM.
Pawnless tablebases with at least one pair of identical pieces of the same
color, such as `KRRNvKRBN`, can be generated on a PC with 64 GB of RAM.
For an 8-piece tablebase with one pawn, 32 GB of RAM should be sufficient.

The generator writes and reads a massive number of temporary files to and
from disk. Use of an SSD is at your own risk (SSDs have limited write
endurance, normally expressed by a TBW rating: *terabytes written*).
Using an HDD should be feasible because the number of disk seeks is expected
to be relatively small (for 7- and 8-piece tables). I suspect that an HDD's
sequential throughput should be sufficient, although this has not yet been
verified.

## Compilation

In the `src` directory, run:

```bash
make tb8gen tb8genp
```

The main dependency is `libzstd` and its header files (`libzstd-devel` on
Fedora).

You may want to set `COMPRESSION_THREADS` to a value no greater than the
number of CPU cores or hardware threads available.

## Usage

To generate `KQvK` (or `KPvK`):

```bash
./tb8gen KQvK
./tb8genp KPvK
```

This creates the files `KQvK.rtbw`, `KQvK.rtbz`, and `KQvK.txt`
(containing statistics) in the current working directory.

During generation, the program can be stopped with `Ctrl+C` (or otherwise
terminated). When restarted later, it will resume where it left off.
(This may not work correctly if the computer is abruptly powered off or
crashes before disk caches are synchronized.)

### Environment Variables

The generator looks in `$RTBPATH` for tablebase files for subtables
(directly or indirectly reachable by capture and/or promotion).
Only WDL files (ending in `.rtbw`) are required.

A temporary working directory with the name of the tablebase (e.g. `KQvK`)
is created in the directory specified by `$TB8DIR`, if set, and otherwise
in the current working directory.

## Command-Line Options

### `--threads n` (or `-t n`)

Use `n` CPU threads to generate and compress the tablebase files.
Currently, this does not affect the number of threads used for compressing
and decompressing temporary data written to and read from disk (this number
is configured in the Makefile).

### `--path /path/to/tbs`

Override the `$RTBPATH` environment variable.
Multiple directories separated by `:` may be specified.

### `--workdir workdirectory`

Specify the location for the temporary working directory
(which will then be `workdirectory/KQvK/`).

### `--layout 0/1/2` (or `-l 0/1/2`)

Specify the layout of the generated tablebase file.

#### Layout 0

Layout 0 corresponds to the layout of existing Syzygy tablebase files and
may result in files compatible with existing chess software. However, the
generator may still choose to use a format extension, rendering the file
incompatible. This option is not available for `tb8genp`.

Layout 0 is the default for tablebases with up to 6 pieces.
(However, Layout 1 may result in smaller files.)

#### Layout 1

Layout 1 is the default for 7-piece tablebase files.
Using this layout, the tablebase is compressed as 10 or 24
(pawnless/pawnful) separate tables
(20 or 48 if the table is two-sided).

#### Layout 2

Layout 2 is the default for 8-piece tablebase files.
This layout corresponds to 462 or 1512 separate compressed tables,
or double that number for two-sided tables.

By compressing the file as multiple tables for 7- and 8-piece tablebases,
the compression algorithm can run entirely in RAM.
For smaller tablebases, the extra overhead of multiple internal tables
increases the total compressed size. For large tablebases, however, this
approach tends to produce smaller files.

### `-c`

Completely clean up and delete the temporary working directory.

### `-g`

Stop after generating the table
(generating bitslices and merging the slices).

## Testing the Tables

Use the `probetool` program (the `new_format` branch) available here:

[https://github.com/syzygy1/probetool/tree/new_format](https://github.com/syzygy1/probetool/tree/new_format)

