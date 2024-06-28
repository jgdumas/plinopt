--------------------------------------------------------------------------------
# PLinOpt: C++ routines handling linear, bilinear & trilinear programs
--------------------------------------------------------------------------------

**Authors**:  Jean-Guillaume Dumas, Bruno Grenet, Clément Pernet, Alexandre Sedoglavic
- [ J-G. Dumas, C. Pernet, A. Sedoglavic; Strassen's algorithm is not optimally accurate, Feb. 2024](https://hal.science/hal-04441653)
- [ J-G. Dumas, B. Grenet; In-place accumulation of fast multiplication formulae, Jul. 2023](https://hal.science/hal-04167499)



**Requirements**:
- C++, pkg-config
- [LinBox](https://linalg.org/), dev: headers & library (LinBox version ≥ 1.7.0; givaro version ≥ 4.2.0)




**Installation**:
- Requires some distribution packages like:
           `sudo apt install git make g++ pkg-config liblinbox-dev`
           (sometimes also `sudo apt install libntl-dev libiml-dev libflint-dev`).
- Then just run `make`, in order to produce the following executable programs
- See also [`bin/auto-docker.run`](https://github.com/jgdumas/plinopt/blob/main/bin/auto-docker.run)
- `make check`, will run correctness programs on the examples files in the `data` directory


**From matrices to programs**:
|  |  |
| :--------- | :------ |
|`bin/optimizer`| produces a small program computing a linear transformation|
|`bin/sparsifier`| factors an MxN matrix into a sparser one, times an NxN matrix |
|`bin/factorizer`| factors an MxN matrix into a sparser MxK, times an KxN matrix |
|`bin/inplacer`| produces an in-place program from a linear transformation|
|`bin/trilplacer`| produces an in-place program from a trilinear transformation|
|  |  |



**Optimizing programs**:
|  |  |
| :--------- | :------ |
|`bin/compacter`| rewrites a simple program using less variables |
|`bin/transpozer`| transposes a program, via Tellegen's transposition principle|
|  |  |



**Tools**:
|  |  |
| :--------- | :------ |
|`bin/matrix-transpose`| transposes a matrix from an SMS file |
|`bin/sms2pretty`| pretty prints a matrix from a file |
|`bin/PMchecker`| asserts correctness of program, with respect to a matrix |
|`bin/MMchecker`| asserts correctness of trilinear program for matrix-multiplication |
|  |  |



**Matrix Syntax**: SMS format (numbering from 1 to m), see [Sparse Integer Matrix Collection](https://hpac.imag.fr)
|  |  |  |  |
| :--- | :--- | :--- | :--- |
| Starts with| `m` | `n` | `'R'` |
| then | `i` | `j` | `value` |
| ends with| `0` | `0` | `0` |
|  |  |  |  |
- `data` directory contains example matrices in this format




**Program Syntax**:
- Line by line
- Operators are: `:=`, `+`, `-`, `*`, `/`, `;`
- Each line not containing `:=` is ignored
- Everything after ';' is ignored
- Everything not an operator, nor spaces, is a `variable`:
	- Constant variables start by `c`
	- Input variables start by `i`
	- Output variables start by `o`
- Program must start by loading each input variable into a temporary `t1:=i1`


**Examples**:
- `./bin/sms2pretty data/Lw.sms data/Rw.sms data/Pw.sms`: pretty print HM representation of Strassen-Winograd's fast 2x2 multiplication algorithm
- `./bin/matrix-transpose data/Pw.sms`: the transposed matrix
- `./bin/MMchecker data/Lw.sms data/Rw.sms data/Pw.sms`: probabilistically checking Strassen-Winograd's matrix-multiplication algorithm
- `./bin/MMchecker data/Lo.sms data/Ro.sms data/Po.sms 32 3 1013`: probabilistically checking matrix-multiplication algorithm using 32 bits random elements and 1013 as a placeholder for sqrt(3)
- `./bin/optimizer data/cyclic.sms`: a program computing that matrix-vector product
- `./bin/optimizer data/Rr.sms -q 3`: a program computing that matrix-vector product modulo 3
- `./bin/transpozer data/test.prg`: a program computing the transposed program
- `./bin/compacter data/test.prg`: a more compact program
- `./bin/optimizer -D data/Pi.sms`: a program computing that matrix-vector product
- `./bin/optimizer data/Li.sms | ./bin/compacter -s`: a compact program computing that matrix-vector product
- `./bin/matrix-transpose data/Pi.sms | ./bin/optimizer -K | ./bin/transpozer`: a program computing that matrix-vector product
- `./bin/sparsifier -c 4 data/Lr.sms`: a factorization of that matrix into a sparser one (also with many 1s) by an alternate change of basis (CoB) 4x4 matrix
- `./bin/factorizer -k 6 data/Lr.sms`: a factorization of that matrix into a 7x6 sparser one by a 6x4 matrix
- `./bin/inplacer data/Lw.sms`: in-place matrix-vector accumulating multiplication
- `./bin/trilplacer data/Lw.sms data/Rw.sms data/Pw.sms`: in-place version of Strassen-Winograd's fast 2x2 accumulating multiplication
- `./bin/trilplacer data/Lk.sms data/Rk.sms data/Pk.sms -e`: in-place version of Karatsuba's fast accumulating polynomial multiplication
- `./bin/optimizer -q 17 data/Lo.sms | ./bin/compacter -s | ./bin/PMchecker -q 17 -M data/Lo.sms`: creates a program for `Lo.sms` modulo 17, compacts it with less variables, then checks consistency
