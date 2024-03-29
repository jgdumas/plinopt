--------------------------------------------------------------------------------
# PLinOpt: a collection of C++ routines handling linear & bilinear programs
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



**About**:
|  |  |
| :--------- | :------ |
|`transpozer`| transposes a program, via Tellegen's transposition principle|
|`optimizer`| produces a small program computing a linear transformation|
|`sparsifier`| factors an MxN matrix into a sparser one, times an NxN matrix |
|`factorizer`| factors an MxN matrix into a sparser MxK, times an KxN matrix |
|`inplacer`| produces an in-place program from a linear transformation|
|`bilplacer`| produces an in-place program from a bilinear transformation|
|  |  |



**Tools**:
|  |  |
| :--------- | :------ |
|`matrix-transpose`| transposes a matrix from an SMS file |
|`sms2pretty`| pretty prints a matrix from an SMS file |
|`MMchecker`| asserts correctness of bilinear program for matrix-multiplication |
|`compacter`| rewrites a simple program using less variables |
|  |  |



**Matrix Syntax**: SMS format, see [Sparse Integer Matrix Collection](https://hpac.imag.fr)
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
- `./sms2pretty data/Lw.sms data/Rw.sms data/Pw.sms`: pretty print HM representation of Strassen-Winograd's fast 2x2 multiplication algorithm`
- `./matrix-transpose data/Pw.sms`: the transposed matrix
- `./MMchecker data/Lw.sms data/Rw.sms data/Pw.sms`: probabilistically checking Strassen-Winograd's matrix-multiplication algorithm
- `./MMchecker data/Lo.sms data/Ro.sms data/Po.sms 32 3 1013`: probabilistically checking matrix-multiplication algorithm using 32 bits random elements and 1013 as a placeholder for sqrt(3)
- `./optimizer data/cyclic.sms`: a program computing that matrix-vector product
- `./optimizer data/Rr.sms -q 3`: a program computing that matrix-vector product modulo 3
- `./transpozer data/test.prg`: a program computing the transposed program
- `./compacter data/test.prg`: a more compact program
- `./optimizer -D data/Pi.sms`: a program computing that matrix-vector product
- `./optimizer data/Li.sms |./compacter -s`: a compact program computing that matrix-vector product
- `./matrix-transpose data/Pi.sms | ./optimizer -K | ./transpozer`: a program computing that matrix-vector product
- `./sparsifier -c 4 data/Lr.sms`: a factorization of that matrix into a sparser one (also with many 1s) by an alternate change of basis(CoB) 4x4 matrix
- `./factorizer -k 6 data/Lr.sms`: a factorization of that matrix into a 7x6 sparser one by a 6x4 matrix
- `./inplacer data/Lw.sms`: in-place matrix-vector accumulating multiplication
- `./bilplacer data/Lw.sms data/Rw.sms data/Pw.sms`: in-place version of Strassen-Winograd's fast 2x2 accumulating multiplication
- `./bilplacer data/Lk.sms data/Rk.sms data/Pk.sms -e`: in-place version of Karatsuba's fast accumulating polynomial multiplication
