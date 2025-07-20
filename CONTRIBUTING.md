# Contributing to MM/GBSA Analysis Package

Thank you for your interest in contributing to the MM/GBSA Analysis Package! This document provides guidelines and information for contributors.

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
git clone https://github.com/your-username/mmgbsa.git
cd mmgbsa

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -e .  # Install in development mode
```

### Development Dependencies
```bash
pip install -r requirements-dev.txt  # If available
# Or install manually:
pip install pytest black flake8 mypy sphinx
```

## 🧪 Testing

### Run Tests
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=mmgbsa

# Run specific test file
pytest test/test_mmgbsa.py

# Run with verbose output
pytest -v
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
mmgbsa/
├── mmgbsa/              # Main package
│   ├── __init__.py
│   ├── core/           # Core functionality
│   ├── analysis/       # Analysis modules
│   ├── utils/          # Utility functions
│   └── config/         # Configuration handling
├── tests/              # Test suite
├── docs/               # Documentation
├── examples/           # Example scripts
├── setup.py           # Package setup
├── requirements.txt   # Dependencies
└── README.md         # Project documentation
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