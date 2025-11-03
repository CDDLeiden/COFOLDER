# Boltz-Tools

## Introduction 
Boltz-Tools is a collection of command-line utilities for performing various co-folding tasks using Boltz. It provides subcommands for predicting single systems, screening ligand/protein libraries, using Boltz as an oracle, and validating system configurations.

## Installation
Install directly from GitHub for newest updates:
```
git clone https://github.com/rlvandenbroek/boltz-tools.git
cd boltz-tools; pip install -e .
```

## Usage
The main command is `boltz-tools`, which supports several subcommands:
```
boltz-tools [-h] [-v] {predict,screen,oracle,validate}
```

### Subcommands
- **predict**: Co-fold a single system using Boltz.  
- **screen**: Co-fold a library using Boltz for virtual screening.  
- **oracle**: Use Boltz as an oracle function for single SMILES predictions.  
- **validate**: Validate Boltz system configuration.

Use the -h flag with any command to see detailed usage:
```
boltz-tools -h
boltz-tools screen -h
```

## Contributing
1. Create a feature branch 
2. Make your changes
3. (Add tests for new functionality -- (applicable when tests are in place))
4. (Ensure all tests pass -- (applicable when tests are in place))
5. Include new features in tutorials      
6. Submit a pull request 

## Package Structure
```
examples/                          # Example input files for tutorials 
legacy/                            # Obsolete files
src/
└── boltz_tools/
    ├── __init__.py
    ├── __main__.py
    ├── cli.py                     # Root CLI 
    ├── helpers/                   # General utilities 
    │   ├── __init__.py
    │   ├── command.py             # Boltz CLI 
    │   ├── conformers.py          # CCD handling
    │   ├── system.py              # System YAML modification
    │   └── utils.py               # General utilities
    ├── tools/                     # Tools / subcommands 
    │   ├── __init__.py
    │   ├── predict.py
    │   ├── validate.py
    │   ├── screen.py
    │   └── oracle.py
    └── analytics/                 # Data analysis utilities 
        ├── __init__.py
        ├── stats.py
        └── plotting.py
tests/                             # Unit tests
tutorials/                         # Tutorials
LISENCE
pyproject.toml
README.md
```
