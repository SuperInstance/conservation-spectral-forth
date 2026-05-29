# Conservation Spectral SDK — Forth

**Stack-Based = Sheaf-Theoretic Spectral Graph Analysis**

The Forth implementation of the Conservation Spectral framework. Builds graph Laplacians, runs power-iteration eigendecomposition, computes conservation ratios, and performs anomaly detection — all through Forth's defining words and stack discipline.

## The Aha Moment

Forth's stack *is* the stalk in sheaf theory. Each word is a **local section**. Composition via concatenation (the only way to combine words in Forth) *is* restriction. `FULL-ANALYSIS` is the **global section** that glues all local sections together. No other language makes this correspondence literal — in Forth, the mathematical structure of sheaves isn't modeled, it *is* the programming model. When you push values onto the stack, transform them, and pop results, you're performing stalk-wise operations. When you compose words, you're building restriction maps. This isn't a metaphor — it's the same algebra.

## How to Use

```bash
# Install Gforth
sudo apt install gforth

# Run the full analysis on the built-in 5-node chord graph
gforth conservation_spectral.fs -e "FULL-ANALYSIS bye"
```

Expected output: eigenvalues, spectral gap, Cheeger constant, conservation ratios, and anomaly detection for the example graph.

### Extending

Define your own graphs by setting `TRANSITION` matrix entries and calling the analysis words:

```forth
5 GRAPH-SIZE !
0E1 0 1 TRANSITION MATRIX-SET   \ edge 0→1
0E1 1 0 TRANSITION MATRIX-SET   \ edge 1→0
\ ... set all edges
FULL-ANALYSIS
```

## Architecture

- **Matrix operations**: `MATRIX-GET`, `MATRIX-SET`, `CLEAR-MATRIX` — row-major float arrays
- **Vector operations**: `VEC-NORM`, `VEC-NORMALIZE`, `VEC-COPY2TO1` — operate on `WORK-VEC`
- **Core analysis**: `BUILD-LAPLACIAN`, `POWER-ITERATION`, `COMPUTE-CONSERVATION`, `DETECT-ANOMALY`
- **Full pipeline**: `FULL-ANALYSIS` — runs everything and prints results
- Max graph size: `MAX-N` (default 10, changeable)

## Connection to the Conservation Spectral Framework

This is one of 12+ implementations of the same spectral conservation analysis across different languages and compute platforms. The framework measures how well a graph's Laplacian eigenvectors conserve an attribute — the **conservation ratio** α(G,a) = (a^T L a) / (λ_max ‖a‖²).

The Forth implementation proves that the entire pipeline — from graph construction through spectral analysis to anomaly detection — requires nothing beyond a stack and dictionary. No objects, no types, no memory management. Just words composing words.

## Related Repos

- [conservation-spectral-lisp](https://github.com/SuperInstance/conservation-spectral-lisp) — Symbolic theorem proving (the other "thinking" language)
- [conservation-spectral-asm](https://github.com/SuperInstance/conservation-spectral-asm) — Register-level AVX2 implementation
- [conservation-spectral-ptx](https://github.com/SuperInstance/conservation-spectral-ptx) — GPU assembly (NVIDIA PTX)
- [conservation-spectral-v2](https://github.com/SuperInstance/conservation-spectral-v2) — Reference Python implementation
- [conservation-geometry](https://github.com/SuperInstance/conservation-geometry) — Geometric visualizations of the framework

## License

MIT
