# Task Completion Checklist

## Before Completing a Task

### Code Quality
- [ ] Code follows project style conventions (119 char line length, 4-space indentation)
- [ ] Code is properly formatted with YAPF: `yapf -i *.py`
- [ ] Code passes linting: `flake8 --max-line-length=119 *.py`
- [ ] Functions and classes have appropriate docstrings
- [ ] Imports are organized correctly (standard, third-party, local)

### Testing
- [ ] New functionality has corresponding tests
- [ ] Tests run successfully: `pytest`
- [ ] Code coverage is maintained: `pytest --cov=mofography`
- [ ] Test files are named appropriately (test_*.py)
- [ ] Tests cover edge cases and error conditions

### Documentation
- [ ] README is updated if user-facing changes were made
- [ ] API documentation is updated for new functions/classes
- [ ] Examples are updated if behavior changed
- [ ] Tutorial notebooks are updated if relevant

### Integration
- [ ] New code integrates properly with existing modules
- [ ] Import statements in __init__.py are updated if needed
- [ ] Dependencies are properly declared if new ones were added
- [ ] Version number is bumped if breaking changes were made

### Performance and Validation
- [ ] Code performs efficiently on typical datasets
- [ ] Memory usage is reasonable for expected data sizes
- [ ] Rust components (if modified) compile and run correctly
- [ ] Integration with Dask for parallel computing works if relevant

## Commands to Run Before Submitting

```bash
# Format code
yapf -i mofography/*.py

# Run linting
flake8 --max-line-length=119 mofography/*.py

# Run tests
pytest
pytest --cov=mofography

# Build documentation (if docs changed)
cd docs && make html

# Check if package imports correctly
python -c "import mofography as mgr; print('Import successful')"
```

## After Task Completion

### Git Operations
- [ ] Stage relevant files: `git add .`
- [ ] Create descriptive commit message
- [ ] Push changes: `git push origin master`
- [ ] Create pull request if working on a branch

### Cleanup
- [ ] Remove temporary files
- [ ] Clean up __pycache__ directories
- [ ] Update CHANGELOG if needed
- [ ] Close any temporary debugging code

### Verification
- [ ] Test the changes work as expected
- [ ] Verify no regressions in existing functionality
- [ ] Check documentation builds correctly
- [ ] Ensure package can be imported and used

## Special Considerations

### For Rust Components
- [ ] Rust code compiles: `cargo build`
- [ ] Rust tests pass: `cargo test`
- [ ] Python-Rust integration works correctly
- [ ] Performance benchmarks are run if relevant

### For Documentation Changes
- [ ] Documentation builds without errors
- [ ] Links are working correctly
- [ ] Images and figures display properly
- [ ] Code examples are tested and work

### For Breaking Changes
- [ ] Deprecation warnings are added if appropriate
- [ ] Migration guide is provided
- [ ] Backward compatibility is considered
- [ ] Major version bump is justified