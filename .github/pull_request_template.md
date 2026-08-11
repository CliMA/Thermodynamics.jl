# Purpose

<!-- What does this PR change, and why? Link any relevant issue. -->

# Content

<!-- A short list of the main changes. -->

-

---

# Checklist

- [ ] Code is formatted (`julia --project=.dev .dev/climaformat.jl .`)
- [ ] Docstrings added or updated for any public function that changed
- [ ] Tests added or updated, and the suite passes (`julia --project=test test/runtests.jl`)
- [ ] `NEWS.md` entry added, with a badge if the change is a bugfix, breaking, or a
      performance change
- [ ] Documentation builds (`julia --project=docs docs/make.jl`), if docs changed

For changes that affect numerical results:

- [ ] The change in behaviour is described in `NEWS.md`
- [ ] Downstream CI (ClimaAtmos, ClimaCoupler, ClimaLand, KinematicDriver) has been
      considered, and any expected differences are noted above
- [ ] Zero-allocation and type-stability tests still pass for the affected kernels
