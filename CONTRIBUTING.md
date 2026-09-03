# Contributing to OpenGBSA

Thank you for your interest in contributing to OpenGBSA! This document provides guidelines and information for contributors.

## 🚀 Quick Start

1. **Fork** the repository
2. **Clone** your fork locally
3. **Create** a feature branch
4. **Make** your changes
5. **Test** your changes
6. **Submit** a pull request

## 📋 Development Setup

### Prerequisites
- Python 3.8+
- Git
- OpenMM 8.0+

### Installation
```bash
# Clone your fork
git clone https://github.com/your-username/opengbsa.git
cd opengbsa

# Create environment (conda recommended -- several dependencies are conda-only)
conda create -n opengbsa python=3.10
conda activate opengbsa
conda install -c conda-forge openmm mdtraj rdkit parmed

# Install the package in development mode
pip install -e ".[dev]"
```

### Development Dependencies
The `dev` extra in `pyproject.toml` installs `pytest`, `black`, `flake8`, `mypy`, and `sphinx`:
```bash
pip install -e ".[dev]"
```

## 🧪 Testing

### Run Tests
```bash
# Run the automated unit test suite
pytest test/unit

# Run with coverage
pytest test/unit --cov=mmgbsa

# Run a specific test file
pytest test/unit/test_config.py

# Run with verbose output
pytest test/unit -v
```

### Test Data
- Test files are in the `test/` directory
- Use small, representative datasets for testing
- Don't commit large trajectory files

## 📝 Code Style

### Python Style Guide
We follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) with some modifications:

- **Line length**: 88 characters (Black default)
- **Docstrings**: Google style
- **Type hints**: Required for public functions

### Code Formatting
```bash
# Format code with Black
black .

# Check code style with flake8
flake8 .

# Type checking with mypy
mypy mmgbsa/
```

### Pre-commit Hooks
```bash
# Install pre-commit hooks
pre-commit install

# Run manually
pre-commit run --all-files
```

## 🔧 Development Workflow

### 1. Issue Reporting
- Use the issue template
- Provide clear, reproducible examples
- Include system information and error messages

### 2. Feature Development
- Create a feature branch: `git checkout -b feature/amazing-feature`
- Write tests for new functionality
- Update documentation
- Ensure all tests pass

### 3. Bug Fixes
- Create a bug fix branch: `git checkout -b fix/bug-description`
- Add regression tests
- Verify the fix works
- Update relevant documentation

### 4. Documentation
- Update README.md if needed
- Add docstrings to new functions
- Update API documentation
- Include usage examples

## 📊 Pull Request Guidelines

### Before Submitting
- [ ] Code follows style guidelines
- [ ] Tests pass locally
- [ ] Documentation is updated
- [ ] No large files are committed
- [ ] Commit messages are clear and descriptive

### Pull Request Template
```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

## Testing
- [ ] Tests pass locally
- [ ] New tests added for new functionality
- [ ] All existing tests still pass

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No breaking changes (or documented if necessary)
```

## 🏗️ Project Structure

```
opengbsa/
├── mmgbsa/              # Main package
│   ├── __init__.py
│   ├── mmgbsa_core.py   # GBSACalculator and core energy pipeline
│   ├── runner.py        # MMGBSARunner: config-driven execution
│   ├── config.py        # ConfigManager: YAML validation/loading
│   ├── topology.py      # Topology loading (Amber/GROMACS/CHARMM)
│   ├── reporting.py     # HTML report generation
│   └── forcefields/     # Bundled force field files
├── test/
│   ├── unit/            # Automated pytest suite (run in CI)
│   ├── configs/         # Self-contained validation datasets
│   └── manual/          # Ad hoc developer scripts, not run in CI
├── docs/                # Sphinx documentation source
├── pyproject.toml       # Package metadata and dependencies
└── README.md            # Project documentation
```

## 🐛 Bug Reports

### Bug Report Template
```markdown
## Bug Description
Clear description of the bug

## Steps to Reproduce
1. Step 1
2. Step 2
3. Step 3

## Expected Behavior
What you expected to happen

## Actual Behavior
What actually happened

## Environment
- OS: [e.g., Ubuntu 20.04]
- Python: [e.g., 3.9.7]
- MM/GBSA Version: [e.g., 0.0.4]
- OpenMM Version: [e.g., 8.0.0]

## Additional Information
Any other relevant information
```

## 💡 Feature Requests

### Feature Request Template
```markdown
## Feature Description
Clear description of the requested feature

## Use Case
Why this feature would be useful

## Proposed Implementation
Optional: How you think it could be implemented

## Alternatives Considered
Optional: Other approaches you've considered
```

## 📚 Documentation

### Building Documentation
```bash
# Install documentation dependencies
pip install sphinx sphinx-rtd-theme

# Build documentation
cd docs
make html

# View documentation
open _build/html/index.html
```

### Documentation Guidelines
- Use clear, concise language
- Include code examples
- Add screenshots for GUI features
- Keep documentation up to date

## 🤝 Community Guidelines

### Code of Conduct
- Be respectful and inclusive
- Help others learn and grow
- Provide constructive feedback
- Follow the project's coding standards

### Communication
- Use GitHub issues for bug reports and feature requests
- Use GitHub discussions for general questions
- Be patient and helpful with new contributors

## 📄 License

By contributing to this project, you agree that your contributions will be licensed under the MIT License.

## 🙏 Acknowledgments

Thank you to all contributors who have helped make this project better!

---

For questions about contributing, please open an issue or contact the maintainers. 