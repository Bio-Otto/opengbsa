# GitHub Repository Setup Guide

## 🚀 Repository Oluşturma

### 1. GitHub'da Yeni Repository
- [GitHub.com](https://github.com) → **Bio-Otto** hesabına gidin
- **"+"** → **"New repository"**
- **Repository name**: `opengbsa`
- **Description**: `OpenGBSA - Molecular Mechanics / Generalized Born Surface Area Analysis Package`
- **Visibility**: Public
- **License**: MIT License
- **Create repository**

### 2. Repository Ayarları

#### Social Preview (Repository Banner)
- **Settings** → **General** → **Social preview**
- **Upload image**: `assets/github_banner.png` (1280x640)
- **Save**

#### Repository Topics
- **About** → **Topics** ekleyin:
  - `molecular-dynamics`
  - `mm-gbsa`
  - `computational-chemistry`
  - `drug-design`
  - `openmm`
  - `mdtraj`
  - `python`
  - `biochemistry`
  - `protein-ligand`
  - `binding-energy`

#### Repository Description
```
OpenGBSA: A comprehensive Molecular Mechanics/Generalized Born Surface Area analysis package with advanced features including entropy analysis, per-residue decomposition, and protein-ligand interaction analysis.
```

### 3. Push Kodu
```bash
git push -u origin main
```

## 📋 Repository Özellikleri

### Badges (README.md'ye eklenecek)
```markdown
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![OpenMM 8.0+](https://img.shields.io/badge/OpenMM-8.0+-green.svg)](https://openmm.org/)
[![Version](https://img.shields.io/badge/version-0.0.1-blue.svg)](https://github.com/Bio-Otto/opengbsa)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXXX-blue.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

### Repository Stats
- **Stars**: 0 (başlangıç)
- **Forks**: 0 (başlangıç)
- **Issues**: 0 (başlangıç)
- **Pull Requests**: 0 (başlangıç)

## 🎯 GitHub Pages

### Dokümantasyon Sitesi
- **Settings** → **Pages**
- **Source**: Deploy from a branch
- **Branch**: main
- **Folder**: /docs
- **Save**

### Özel Domain (Opsiyonel)
- **Custom domain**: `opengbsa.bio-otto.com`
- **Enforce HTTPS**: ✓

## 📊 Insights

### Traffic Analytics
- **Settings** → **General** → **Features**
- **Issues**: ✓ Enable
- **Pull requests**: ✓ Enable
- **Wikis**: ✓ Enable
- **Discussions**: ✓ Enable

### Security
- **Settings** → **Security**
- **Dependency graph**: ✓ Enable
- **Dependabot alerts**: ✓ Enable
- **Code scanning**: ✓ Enable

## 🔧 GitHub Actions

### CI/CD Pipeline
`.github/workflows/ci.yml` dosyası oluşturun:

```yaml
name: CI

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, "3.10", "3.11"]

    steps:
    - uses: actions/checkout@v3
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v3
      with:
        python-version: ${{ matrix.python-version }}
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install pytest
    - name: Run tests
      run: |
        pytest test/
```

## 📈 Release Management

### v0.0.1 Release
- **Releases** → **Create a new release**
- **Tag version**: `v0.0.1`
- **Release title**: `OpenGBSA v0.0.1 - Initial Release`
- **Description**:
```markdown
## 🎉 Initial Release

### Features
- Complete MM/GBSA energy calculation
- Multiple GB models (OBC2, OBC1, HCT, GBn, GBn2)
- Normal Mode Analysis entropy calculations
- Per-residue energy decomposition
- Protein-ligand interaction analysis (ProLIF)
- YAML-based configuration system
- Advanced visualization and reporting
- Professional logo and branding

### Installation
```bash
pip install opengbsa
```

### Quick Start
```bash
opengbsa --create-config
opengbsa config.yaml
```
```

## 🌟 Community

### Contributing Guidelines
- **CONTRIBUTING.md** dosyası zaten mevcut
- **Code of Conduct** ekleyin
- **Issue templates** oluşturun

### Issue Templates
`.github/ISSUE_TEMPLATE/bug_report.md`:
```markdown
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run command '...'
2. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Environment:**
 - OS: [e.g. Ubuntu 20.04]
 - Python version: [e.g. 3.9]
 - OpenGBSA version: [e.g. 0.0.1]
```

## 📞 Support

### Contact Information
- **GitHub Issues**: [opengbsa/issues](https://github.com/Bio-Otto/opengbsa/issues)
- **Email**: bio.otto@marmara.edu.tr
- **Website**: [compbio-hub.com](https://compbio-hub.com)

### Documentation
- **README.md**: Ana dokümantasyon
- **Wiki**: Detaylı kullanım kılavuzu
- **Examples**: Örnek konfigürasyonlar
- **API Reference**: Kod dokümantasyonu

## 🎯 Sonraki Adımlar

1. **Repository oluşturun** (yukarıdaki adımları takip edin)
2. **Kodu push edin**: `git push -u origin main`
3. **Release oluşturun**: v0.0.1
4. **Dokümantasyon ekleyin**: Wiki ve examples
5. **Community oluşturun**: Issues, discussions, contributions
6. **Promote edin**: Twitter, LinkedIn, academic networks

---

**Not**: Bu rehber repository oluşturulduktan sonra güncellenebilir. 