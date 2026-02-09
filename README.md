# ExtendedStim: A Python Package Addressing both Fermionic and Bosonic Quantum Error Correction Simultaneously


[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.12%2B-blue)](https://www.python.org/)
[![Documentation](https://img.shields.io/badge/docs-available-brightgreen)](docs/)

This program is developed with Python 3.12+ and is mainly used for constructing and testing quantum error-correcting codes and quantum circuits.

## 🔨 1 Project Dependencies

Required:

- [QuTiP](https://qutip.org/) - Quantum toolbox
- [Stim](https://github.com/quantumlib/Stim) - Quantum error-correction simulator
- [Stimbposd](https://github.com/quantumlib/Stim/blob/main/docs/bposd.md) - Stim-based BPOSD decoder
- [Galois](https://galois.readthedocs.io/) - Algebra over $\mathbb{F}_2$
- [NumPy](https://numpy.org/) - Numerical computing library
- [SciPy](https://scipy.org/) - Scientific computing library
- [Matplotlib](https://matplotlib.org/) - Plotting library
- [Qiskit](https://qiskit.org/) - Quantum circuit drawing
- [Mip](https://www.mipengine.org/) - Integer programming solver for computing code distance

Optional:

[tesseract-decoder](https://github.com/quantumlib/tesseract-decoder) - Tesseract decoder; if not installed, Stimbposd is used instead.

## 📖 2 Basic Workflow

### 2.1 Compute code parameters

1. Construct a quantum error-correcting code
2. Compute the code parameters of the quantum error-correcting code

### 2.2 Compute logical error rate

1. Construct a quantum circuit
2. Run Monte Carlo simulations and compare predicted correctness to obtain the logical error rate

## 📄 3 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 4 Contact

- **Author**: Moke
- **Email**: Moke2001@whu.edu.cn
- **Address**: Room S219, Meng Minwei Science and Technology Building, Tsinghua University, Haidian District, Beijing
- **Phone**: +86 130-3373-6868
