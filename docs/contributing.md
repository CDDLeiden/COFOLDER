# Contributing to COFOLDER

Thank you for your interest in contributing to COFOLDER! This guide will help you get started.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/CDDLeiden/COFOLDER.git
   cd COFOLDER
   ```
3. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Setup

Install the supported test, tutorial, documentation, and development tools:

```bash
python -m pip install -e ".[analysis,docs,test,tutorials,development]"
```

If your work changes runner integrations, backend packaging, or the shared CLI workflows, use the backend-specific clean-install acceptance lane described in [Backend Acceptance Tutorial](tutorials/backend-acceptance.md). That lane is intentionally manual and expensive, so it complements routine tests instead of replacing them.

## Development Workflow

### 1. Make Your Changes

- Write clear, documented code
- Follow existing code style and conventions
- Add docstrings for new functions and classes
- Use type hints where appropriate

### 2. Add Tests

When tests are in place:
- Write tests for new functionality
- Ensure existing tests still pass
- Aim for good test coverage

The supported suite is divided into four exhaustive, non-overlapping lanes:

```bash
python scripts/run_test_lane.py core
python scripts/run_test_lane.py contracts-tutorial
python scripts/run_test_lane.py artifact
python scripts/run_test_lane.py acceptance
```

Use `python scripts/run_test_lane.py all` to run their union. The lane runner fails
before pytest starts if a supported `test_*.py` module has missing or overlapping
ownership. Install `.[analysis,test]` for core tests, `.[test,tutorials]` for the combined
contracts/tutorial lane, `.[test,development]` for artifact checks, and
`.[test,acceptance]` for the lightweight acceptance-contract tests.

The artifact lane's end-to-end build and installed-package smoke check is:

```bash
python scripts/verify_release_artifacts.py
```

It builds an sdist and then builds the wheel from the unpacked sdist. It verifies
the supported inventory and UI exclusion before installing the wheel outside the
checkout and exercising metadata, entry-point help, and example copying.

For backend-, runner-, or CLI-adjacent changes:
- keep routine targeted tests in place
- run the backend acceptance notebooks before promoting the change toward `main`
- use one fresh environment per backend (`boltz1`, `boltz2`, `boltz-community`, or `openfold3`) instead of reusing a mixed development environment

Routine CI does not install those backend extras, download models, use external MSA
services, require a GPU, or perform inference. Its `acceptance` lane contains only
lightweight mocked contract checks. Actual Boltz1, Boltz2, Boltz Community, and
OpenFold3 acceptance remains manual and each backend must use a separate environment;
a routine CI pass is not evidence that any backend acceptance run passed.

### 3. Update Documentation

- Update relevant documentation in `docs/`
- Add docstrings following NumPy style
- Update examples if needed
- Add tutorials for new features

If you are integrating a new runner or changing runner behavior, update both:
- [Adding New Runners](tutorials/runners.md)
- [Backend Acceptance Tutorial](tutorials/backend-acceptance.md)

### 4. Submit Pull Request

1. Commit your changes with clear messages:
   ```bash
   git add .
   git commit -m "Add feature: description"
   ```

2. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

3. Create a Pull Request on GitHub targeting the `dev` branch

Before changes that affect runners, backend packaging, or shared CLI behavior are promoted toward `main`, run the backend acceptance lane from [Backend Acceptance Tutorial](tutorials/backend-acceptance.md). These checks are not CI-default and should be treated as a deliberate pre-promotion safeguard.

## Code Style

Run the repository-owned lint baseline from the repository root:

```bash
ruff check src tests scripts tutorials examples
```

### Python Style

- Follow PEP 8 guidelines
- Use meaningful variable names
- Keep functions focused and concise
- Add type hints for function signatures

### Docstring Style

Use NumPy-style docstrings:

```python
def my_function(param1: str, param2: int = 5) -> bool:
    """
    Short description of the function.

    Longer description if needed, explaining the function's
    purpose and behavior in detail.

    Parameters
    ----------
    param1 : str
        Description of param1.
    param2 : int, optional
        Description of param2 (default is 5).

    Returns
    -------
    bool
        Description of return value.

    Examples
    --------
    >>> my_function("test", 10)
    True
    """
    pass
```

## Documentation

Build documentation locally to preview changes:

```bash
mkdocs serve
```

Then visit `http://127.0.0.1:8000` in your browser.

## Pull Request Guidelines

Good pull requests include:

- Clear description of changes
- Reference to related issues
- Updated tests (when applicable)
- Updated documentation
- Clean commit history

### PR Checklist

- [ ] Code follows project style
- [ ] Docstrings added/updated
- [ ] Tests added/updated (when applicable)
- [ ] Documentation updated
- [ ] Tutorials updated for new features
- [ ] No merge conflicts with `dev` branch

## Reporting Issues

When reporting issues, include:

- COFOLDER version
- Python version
- Operating system
- Complete error message
- Minimal example to reproduce
- Expected vs actual behavior

## Feature Requests

Feature requests are welcome! Please:

- Check existing issues first
- Clearly describe the use case
- Explain why it would be useful
- Provide examples if possible

## Questions?

If you have questions about contributing:

- Check existing documentation
- Look at similar PRs
- Open a discussion on GitHub

## Code of Conduct

Be respectful and constructive in all interactions. We aim to maintain a welcoming and inclusive community.

## License

By contributing, you agree that your contributions will be licensed under the same license as the project.
