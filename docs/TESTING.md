# VCF Analysis Agent - Testing Documentation

## Overview

The VCF Analysis Agent maintains a comprehensive, multi-layered testing strategy to ensure reliability, correctness, and maintainability. Our testing approach covers unit testing, integration testing, and end-to-end (E2E) testing with industry-leading coverage standards.

## Testing Philosophy

### Quality Standards
- **86% Unit Test Coverage**: Industry-leading coverage across all core modules
- **100% Critical Path Coverage**: All essential workflows thoroughly tested
- **Comprehensive Error Handling**: Edge cases and failure scenarios covered
- **Real-World Validation**: Tests use actual VCF files and realistic scenarios

### Testing Pyramid
```
    /\     E2E Tests (45 tests)
   /  \    CLI Interface & User Workflows
  /____\   
 /      \  Integration Tests (116 tests)
/        \ Component Interactions & Workflows
\________/
 \      /  Unit Tests (86% coverage)
  \____/   Individual Modules & Functions
```

## Test Suite Structure

### 1. Unit Tests (`tests/unit/`)

**Purpose**: Test individual modules and functions in isolation with extensive mocking.

#### Coverage by Module:
- **Agent Module**: 77% coverage, 38 test cases
  - Tool function testing
  - Model initialization
  - Agent creation and configuration
  
- **CLI Module**: 94% coverage
  - Command parsing and validation
  - Argument handling
  - Error scenarios
  
- **LanceDB Integration**: 94% coverage
  - Database operations
  - Embedding management
  - Query functionality
  
- **GATK Integration**: 100% coverage
  - Validation workflows
  - Tool integration
  - Error handling
  
- **Metrics Module**: 77% coverage, 33 test cases
  - Prometheus metrics
  - Performance tracking
  - Error monitoring
  
- **Graph Integration**: 88% coverage
  - Kuzu database operations
  - Graph queries
  - Data modeling

#### Key Test Files:
- `test_agent_comprehensive.py`: Core agent functionality
- `test_cli_comprehensive.py`: CLI interface testing
- `test_lancedb_integration_comprehensive.py`: Vector database testing
- `test_metrics_comprehensive.py`: Monitoring and metrics
- `test_gatk_integration_unit.py`: GATK tool integration
- `test_graph_integration_comprehensive.py`: Graph database testing

### 2. Integration Tests (`tests/integration/`)

**Purpose**: Test component interactions and end-to-end workflows.

#### VCF Processing Workflows (29 tests)
- **Validation Workflows**: VCF → validation → results
- **BCFtools Tool Chains**: view → query → stats workflows
- **Agent Tool Workflows**: Tool chain execution and graph database integration
- **File I/O Workflows**: Temporary file handling and multi-file processing
- **Error Handling**: Cascading errors, partial failures, recovery
- **Performance**: Large file handling, concurrent processing simulation
- **External Systems**: Database and file system integration

#### VCF Comparison Workflows (11 tests)
- **Basic Comparisons**: Comparison workflows with preprocessing
- **Normalization**: Normalization and comparison workflows
- **Multi-Sample**: Multiple file and batch validation workflows
- **Temporary Files**: Compatible VCF creation and comparison

#### Additional Integration Tests (76 tests)
- CLI agent tracing and metrics
- API endpoint testing
- Database integration scenarios
- Error propagation testing

### 3. End-to-End (E2E) Tests (`tests/integration/e2e/`)

**Purpose**: Validate complete user workflows and CLI interface functionality.

#### CLI VCF Processing Tests (29 tests)

**VCF Validation Tests (4 tests)**:
- Valid VCF file validation via CLI
- Invalid VCF file handling
- Non-existent file error handling
- Multiple VCF file validation

**BCFtools Operations Tests (6 tests)**:
- Header viewing through agent interface
- Statistics generation via CLI
- Variant querying with format specification
- Filtering operations with quality thresholds
- Normalization with reference genome
- Error handling for invalid files

**VCF Comparison Tests (3 tests)**:
- Basic VCF comparison workflows
- Self-comparison validation
- Error handling for missing files

**VCF Analysis Tests (3 tests)**:
- Summary analysis generation
- Detailed analysis workflows
- Graph database integration

**Workflow Integration Tests (3 tests)**:
- Validate-then-analyze workflows
- Stats-then-summary workflows
- Multi-step analysis pipelines

**Error Handling Tests (5 tests)**:
- Empty command handling
- Invalid command processing
- Malformed file path handling
- Very long command processing
- Special character handling

**Performance Tests (3 tests)**:
- File handling capabilities
- Concurrent usage simulation
- Timeout management

**Temporary File Tests (2 tests)**:
- Temporary VCF file operations
- Temporary directory workflows

#### CLI Direct Commands Tests (16 tests)

**Direct Command Tests (5 tests)**:
- Kuzu database population
- Missing required arguments validation
- Custom database path handling
- Non-existent file error handling
- Invalid VCF file processing

**Argument Parsing Tests (7 tests)**:
- Help command functionality
- No arguments handling
- Invalid subcommand processing
- Model selection arguments
- Raw mode flag handling
- Credentials argument processing
- Reference FASTA argument handling

**LanceDB Command Tests (3 tests)**:
- Database initialization with custom paths
- Variant addition with minimal arguments
- Embedding search functionality

**Error Handling Tests (1 test)**:
- Invalid JSON updates handling
- Invalid embedding format processing
- Missing database path handling
- Invalid SQL filter processing

## Running Tests

### Basic Test Execution

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run specific test categories
pytest tests/unit/          # Unit tests only
pytest tests/integration/   # Integration tests only
pytest tests/integration/e2e/  # E2E tests only
```

### Coverage Testing

```bash
# Run tests with coverage
pytest --cov=src/vcf_agent --cov-report=term-missing

# Generate HTML coverage report
pytest --cov=src/vcf_agent --cov-report=html

# Coverage with specific threshold
pytest --cov=src/vcf_agent --cov-fail-under=85
```

### Docker Testing

#### Container Testing

```bash
# Build and test in container
docker-compose --profile development up -d
docker-compose exec vcf-agent-dev pytest

# Run specific test categories in container
docker-compose exec vcf-agent-dev pytest tests/unit/ -v
docker-compose exec vcf-agent-dev pytest tests/integration/ -v

# Run tests with coverage in container
docker-compose exec vcf-agent-dev pytest --cov=src/vcf_agent --cov-report=term-missing
```

#### Multi-Architecture Testing

```bash
# Test AMD64 build
./scripts/docker-build.sh --platform linux/amd64 --target development
docker run --rm vcf-analysis-agent:dev-amd64 pytest

# Test ARM64 build (if on compatible hardware)
./scripts/docker-build.sh --platform linux/arm64 --target development
docker run --rm vcf-analysis-agent:dev-arm64 pytest
```

#### Container Integration Testing

```bash
# Test full stack with observability
docker-compose up -d
docker-compose exec vcf-agent python -m vcf_agent.cli --help

# Test VCF ingestion in container
docker-compose exec vcf-agent python -m vcf_agent.cli ingest-vcf --vcf-file sample_data/minimal.vcf.gz

# Test AI analysis in container
docker-compose exec vcf-agent python -m vcf_agent.cli ask "vcf_analysis_summary_tool: sample_data/minimal.vcf.gz"

# Verify observability stack
curl -f http://localhost:9090/-/healthy  # Prometheus
curl -f http://localhost:3000/api/health  # Grafana
curl -f http://localhost:16686/  # Jaeger
```

#### Performance Testing in Containers

```bash
# Test resource limits
docker run --memory="4g" --cpus="2" vcf-analysis-agent:latest pytest tests/integration/

# Test with large files (if available)
docker run -v $(pwd)/large_data:/app/large_data vcf-analysis-agent:latest \
  python -m vcf_agent.cli ingest-vcf --vcf-file /app/large_data/large.vcf.gz

# Monitor resource usage during tests
docker stats vcf_analysis_agent
```

#### Security Testing

```bash
# Test non-root execution
docker run --rm vcf-analysis-agent:latest id
# Should show: uid=10001(appuser) gid=10001(appuser)

# Test file permissions
docker run --rm vcf-analysis-agent:latest ls -la /app/

# Vulnerability scanning
./scripts/docker-build.sh --scan
trivy image vcf-analysis-agent:latest
```

### Selective Test Execution

```bash
# Run tests excluding Kuzu (due to segfaults)
pytest -k "not kuzu"

# Run specific test files
pytest tests/unit/test_agent_comprehensive.py
pytest tests/integration/e2e/test_cli_vcf_processing.py

# Run tests matching pattern
pytest -k "test_cli"
pytest -k "validation"
```

### Performance and Debugging

```bash
# Run tests with timing information
pytest --durations=10

# Run tests with detailed output
pytest -v -s

# Run tests with debugging
pytest --pdb

# Run tests in parallel (if pytest-xdist installed)
pytest -n auto
```

## Test Infrastructure

### CLI Test Runner

The E2E CLI tests use a sophisticated test runner that:
- Executes CLI commands via subprocess
- Manages environment variables and PYTHONPATH
- Implements timeout handling (60 seconds per operation)
- Provides proper error handling and recovery
- Supports both agent "ask" commands and direct CLI subcommands

### Test Data Management

#### VCF Test Files
- `sample_test_data/small_valid.vcf`: Minimal valid VCF for basic testing
- `sample_test_data/sample1.vcf`, `sample2.vcf`: Multi-sample testing
- `sample_test_data/detailed.vcf`: More complex VCF structures
- `sample_test_data/22.fa`: Reference genome for normalization testing

#### Temporary File Handling
- Automatic cleanup of temporary files and directories
- Proper isolation between test runs
- Support for concurrent test execution

### Mocking Strategy

#### External Dependencies
- **Kuzu Database**: Mocked for unit tests, real for integration
- **LanceDB**: Mocked for unit tests, temporary instances for integration
- **OpenTelemetry**: Mocked to avoid external dependencies
- **File System**: Temporary directories and files for isolation

#### AI/LLM Integration
- Mock responses for predictable testing
- Real LLM integration for specific integration tests
- Configurable model selection for testing different providers

## Test Quality Assurance

### Coverage Standards
- **Minimum Coverage**: 80% for all modules
- **Current Achievement**: 86% overall coverage
- **Critical Path Coverage**: 100% for essential workflows
- **Error Path Coverage**: Comprehensive edge case testing

### Test Reliability
- **Deterministic Results**: Tests produce consistent outcomes
- **Isolation**: Tests don't interfere with each other
- **Cleanup**: Proper resource cleanup after test execution
- **Timeout Handling**: Prevents hanging tests

### Performance Considerations
- **Test Execution Time**: Optimized for developer productivity
- **Resource Usage**: Efficient memory and CPU utilization
- **Parallel Execution**: Support for concurrent test runs
- **CI/CD Optimization**: Fast feedback in continuous integration

## Continuous Integration

### GitHub Actions Integration
```yaml
# Example CI configuration
- name: Run Tests
  run: |
    pytest --cov=src/vcf_agent --cov-report=xml
    
- name: Upload Coverage
  uses: codecov/codecov-action@v3
  with:
    file: ./coverage.xml
```

### Coverage Reporting
- **HTML Reports**: Generated in `htmlcov/` directory
- **XML Reports**: Generated as `coverage.xml` for CI/CD
- **Terminal Reports**: Immediate feedback during development

### Quality Gates
- **Minimum Coverage**: 80% required for PR approval
- **Test Passing**: All tests must pass before merge
- **Performance**: Test execution time monitoring

## Troubleshooting

### Common Issues

#### Kuzu Segmentation Faults
```bash
# Skip Kuzu tests if experiencing segfaults
pytest -k "not kuzu"
```

#### CLI Timeout Issues
```bash
# Increase timeout for slow operations
export CLI_TEST_TIMEOUT=120
pytest tests/integration/e2e/
```

#### Missing Test Data
```bash
# Ensure test data files exist
ls sample_test_data/
# Download missing files if needed
```

#### Environment Issues
```bash
# Reset test environment
rm -rf .pytest_cache/
rm -rf htmlcov/
rm coverage.xml
```

### Debugging Test Failures

#### Verbose Output
```bash
pytest -v -s tests/path/to/failing_test.py
```

#### Debug Mode
```bash
pytest --pdb tests/path/to/failing_test.py
```

#### Coverage Analysis
```bash
pytest --cov=src/vcf_agent --cov-report=term-missing
```

## Best Practices

### Writing Tests

#### Unit Tests
- Test one function/method per test
- Use descriptive test names
- Mock external dependencies
- Test both success and failure cases
- Maintain test isolation

#### Integration Tests
- Test component interactions
- Use real data when possible
- Test error propagation
- Validate end-to-end workflows
- Clean up resources

#### E2E Tests
- Test from user perspective
- Use realistic scenarios
- Test complete workflows
- Validate CLI interface
- Handle timeouts gracefully

### Test Maintenance

#### Regular Updates
- Update tests when code changes
- Maintain test data currency
- Review coverage reports
- Update documentation

#### Performance Monitoring
- Monitor test execution time
- Optimize slow tests
- Maintain CI/CD efficiency
- Profile resource usage

## Future Enhancements

### Planned Improvements
- **Golden File Testing**: Reference output validation
- **SAMspec Compliance**: Specification adherence testing
- **Performance Benchmarking**: Baseline performance metrics
- **Mutation Testing**: Code quality validation
- **Property-Based Testing**: Automated test case generation

### Infrastructure Enhancements
- **Test Parallelization**: Faster test execution
- **Test Reporting**: Enhanced reporting and analytics
- **Test Data Management**: Automated test data generation
- **CI/CD Integration**: Advanced pipeline integration

---

*This documentation is maintained alongside the test suite and updated with each major testing milestone.* 