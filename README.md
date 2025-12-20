# Boltz-Lab

## Introduction
Boltz-Lab is a collection of command-line utilities for performing various co-folding tasks using Boltz. It provides subcommands for predicting single systems, screening ligand/protein libraries, using Boltz as an oracle, and evaluating system configurations.

## Installation
Install directly from GitHub for newest updates:
```
git clone https://github.com/CDDLeiden/boltz-lab.git
cd boltz-lab; pip install -e .
```

## Usage
The main command is `boltz-lab`, which supports several subcommands:
```
boltz-lab [-h] [-v] {predict,screen,oracle,evaluate}
```

### Subcommands
- **predict**: Co-fold a single system using Boltz.
- **screen**: Co-fold a library using Boltz for virtual screening.
- **oracle**: Use Boltz as an oracle function for single SMILES predictions.
- **evaluate**: Evaluate Boltz system configuration.

Use the -h flag with any command to see detailed usage:
```
boltz-lab -h
boltz-lab screen -h
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
examples/
src/
└── boltz_lab/
    ├── __init__.py
    ├── __main__.py
    ├── cli.py
    ├── recipes/
    │   ├── __init__.py
    │   ├── predict.py
    │   ├── evaluate.py
    │   ├── screen.py
    │   └── oracle.py
    └── modules/
        ├── input/
        │   ├── __init__.py
        │   ├── command.py
        │   └── system.py
        ├── entities/
        │   ├── __init__.py
        │   └── ligand.py
        ├── utils/
        │   ├── __init__.py
        │   ├── helpers.py
        │   └── log.py
        ├── analytics/
        │   └── __init__.py
        └── runners/
            └── __init__.py
docs/
tests/
tutorials/
LICENSE
pyproject.toml
README.md
```
