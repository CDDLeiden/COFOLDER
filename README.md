# COFOLDER

## Introduction
COFOLDER is a collection of command-line utilities for performing various co-folding tasks using Boltz. It provides subcommands for validating single systems, screening ligand/protein libraries, and using Boltz as an oracle.


## Installation
Install directly from GitHub for newest updates:
```
git clone https://github.com/CDDLeiden/cofolder.git
cd cofolder; pip install -e .
```

## Usage
The main command is `cofolder`, which supports several subcommands:
```
cofolder [-h] [-v] {validate,screen,oracle}
```

### Subcommands
- **validate**: Co-fold and validate a single system using Boltz.
- **screen**: Co-fold a library using Boltz for virtual screening.
- **oracle**: Use Boltz as an oracle function for single SMILES predictions.

Use the -h flag with any command to see detailed usage:
```
cofolder -h
cofolder screen -h
```

## Package Structure
```
examples/
src/
└── cofolder/
    ├── __init__.py
    ├── __main__.py
    ├── cli.py
    ├── recipes/
    │   ├── __init__.py
    │   ├── validate.py
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
tests/
tutorials/
LISENCE
pyproject.toml
README.md
```

