# Contributing to OpenGBSA

Bug reports, questions and pull requests are welcome.

* **Report a problem:** open an issue with the command, the configuration file
  and the log (`opengbsa config.yaml 2>&1 | tee run.log`).
* **Development setup:** install from source in editable mode and run the
  tests, as described in the documentation:
  <https://opengbsa.readthedocs.io/en/latest/contributing.html>.

```bash
git clone https://github.com/Bio-Otto/opengbsa.git
cd opengbsa
pip install --no-deps -e .        # in an environment with the dependencies
python -m pytest test/unit -q
```

* Add a regression test under `test/unit/` for every fix and feature; tests
  must not download data.
* Describe the change in `CHANGELOG.md` under "Unreleased".
* A change that alters computed energies needs evidence on the validation
  cases: `python test/configs/check_validation.py <case>`.
