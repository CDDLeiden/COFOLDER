# Contributing to Boltz-Lab

Thank you for your interest in contributing to Boltz-Lab! This guide will help you get started.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/CDDLeiden/boltz-lab.git
   cd boltz-lab
   ```
3. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Setup

Install in development mode with all dependencies:

```bash
pip install -e ".[docs,test]"
```

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

### 3. Update Documentation

- Update relevant documentation in `docs/`
- Add docstrings following NumPy style
- Update examples if needed
- Add tutorials for new features

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

## Code Style

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

- Boltz-Lab version
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
