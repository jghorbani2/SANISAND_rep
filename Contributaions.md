# Contributing to SANISAND repository

Thank you for your interest in contributing to the SANISAND Model repository! This document provides guidelines and instructions for contributing to this project.

## Code of Conduct

Please be respectful and considerate of others when contributing to this project. We aim to foster an inclusive and welcoming environment for everyone.

## How to Contribute

There are several ways you can contribute to this project:

1. Reporting bugs
2. Suggesting enhancements, this includes implementation of more advanced version of SANISAND
3. Submitting code improvements
4. Improving documentation
5. Sharing validation examples

## Reporting Bugs

If you find a bug in the implementation, please open an issue on GitHub with the following information:

- A clear and descriptive title
- A detailed description of the bug
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Any error messages or logs
- Your environment (OS, compiler version, Eigen version)

## Suggesting Enhancements

If you have ideas for enhancements or new features, please open an issue on GitHub with:

- A clear and descriptive title
- A detailed description of the enhancement
- The motivation for the enhancement
- Any potential implementation ideas

## Pull Request Process

1. Fork the repository
2. Create a new branch for your feature or bug fix
3. Make your changes
4. Add or update tests as necessary
5. Ensure all tests pass
6. Update documentation if necessary
7. Submit a pull request

### Code Style Guidelines

Please follow these style guidelines when contributing code:

- Follow C++ best practices
- Document your code with clear comments
- Use meaningful variable and function names
- Write clean, readable code

### Commit Messages

Write clear and meaningful commit messages that explain what changes were made and why. Follow this format:

```
Short (50 chars or less) summary of changes

More detailed explanatory text, if necessary. Wrap it to about 72
characters. The blank line separating the summary from the body is
critical.

Further paragraphs come after blank lines.

- Bullet points are okay too
- Use a hyphen or asterisk followed by a space
```

## Testing

Before submitting a pull request, please ensure that:

1. All existing tests pass
2. You've added tests for new functionality
3. The code compiles without warnings

To run tests:

```bash
cd build
cmake ..
make
ctest
```

## Documentation

If you're adding new features or changing existing ones, please update the relevant documentation. This includes:

- Code comments
- README.md (if applicable)
- User manual
- Technical documentation

## Development Environment Setup

To set up a development environment:

1. Clone the repository
2. Install dependencies:
   - C++ compiler (GCC, Clang, or MSVC)
   - Eigen 3.3+
   - CMake 3.10+
3. Build the project:
   ```bash
   mkdir build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Debug ..
   make
   ```

## Adding New Features

When adding new features, consider:

1. Will this feature be useful to many users?
2. Is it consistent with the overall design?
3. Has it been thoroughly tested?
4. Does it follow best practices?

## Versioning

This project follows semantic versioning (SemVer). Version numbers are in the format MAJOR.MINOR.PATCH:

- MAJOR: Incompatible API changes
- MINOR: Added functionality (backwards compatible)
- PATCH: Bug fixes (backwards compatible)

## Contact

If you have questions or need help, please:

1. Open an issue on GitHub
2. Contact the maintainers directly at javad.ghorbani@gmail.com

Thank you for contributing to the SANISAND Model implementation!
