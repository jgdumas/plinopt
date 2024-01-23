--------------------------------------------------------------------------------
# PLinOpt: a collection of C++ routines handling linear programs
--------------------------------------------------------------------------------

**Authors**:  Jean-Guillaume Dumas


**Requirements**:
- C++
- [LinBox](https://linalg.org/), dev: headers & library



**Automatic linux install & run first benchmarks**:
- Fetch and run [auto-vm.run](https://raw.githubusercontent.com/jgdumas/vespo/main/auto-vm.run)

	- Requires a linux virtual machine or sudoer rights to install packages.
	- Will install distribution packages: `g++`, `make`, `liblinbox-dev`.
	- Then clone and install the `libvespo` libraries.
	- Then run one small example.



**About**:
- `transpozer`: transposes a program, via Tellegen's transposition principle
- `optimizer` : produces a small program computing a linear transformation
- `inplacer`  : produces an in-place program from a bilinear transformation



**Matrix Syntax**:
- SMS format, see [Sparse Integer Matrix Collection](https://hpac.imag.fr)

	- Starts with: `m n 'R'`
	- then: `i j value`
	- ends with: `0 0 0`



**Program Syntax**:
- Line by line
- Operators are: `:=`, `+`, `-`, `*`, `/`, `;`
- Each line not containing `:=` is ignored
- Everything after ';' is ignored
- Everything not an operator, nor spaces, is a `variable`:
	- Constant variables start by `c`
	- Input variables start by `i`
	- Output variables start by `o`
