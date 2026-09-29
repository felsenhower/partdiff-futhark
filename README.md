# partdiff-futhark

A port of [partdiff](https://github.com/parcio/partdiff) that uses [Futhark](https://futhark-lang.org/) for the calculation-heavy stuff.

## Usage

To build this project, you need the [Futhark compiler](https://futhark.readthedocs.io/en/stable/installation.html).

```shell
$ git clone https://github.com/felsenhower/partdiff-futhark
$ cd partdiff-futhark
$ make
$ ./partdiff 1 2 100 1 2 100
```

By passing the `FUTHARK_BACKEND` variable to `make`, you can customize how the futhark library is built.
E.g. by running `make FUTHARK_BACKEND=opencl`, you can built `partdiff_futhark` with OpenCL.
Some options are
- `multicore`: thread-parallel execution, the default
- `c`: sequential C
- `cuda`: CUDA
- `opencl`: OpenCL
- `hip`: HIP
- ...

Run `futhark --help` to get an overview of the back-ends.

## Correctness

This project uses [partdiff_tester](https://github.com/parcio/partdiff_tester) via CI to ensure that the output matches the reference implementation.
It passes the correctness tests with `--strictness=4` (exact match) and with `--valgrind` (i.e. it has no memory leaks).

## Performance

TODO...