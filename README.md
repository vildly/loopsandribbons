# LoopsAndRibbons API 🧬

A RESTful API service for predicting missing loops in protein structures using both simple interpolation and Modeller-based approaches.

## Overview 🔬

This API service provides endpoints for identifying and predicting missing regions (loops) in protein structures. It offers two prediction methods:
1. **Simple Predictor**: Uses linear interpolation with random variation for quick predictions
2. **Modeller Predictor**: Uses Modeller for more accurate, physics-based predictions

## Installation 🛠️

```bash
# Install dependencies using uv
uv pip install -e .
```

### Requirements
- Python 3.8+
- FastAPI
- Biopython
- NumPy
- Modeller (optional, for Modeller-based predictions)

## API Endpoints 🔌

### Basic Usage

```python
import requests

# Initialize a prediction job
response = requests.post(
    "http://localhost:8000/predict",
    json={
        "pdb_file": "path/to/structure.pdb",
        "method": "simple",  # or "modeller"
        "num_conformations": 5
    }
)

# Get prediction results
job_id = response.json()["job_id"]
results = requests.get(f"http://localhost:8000/results/{job_id}")
```

### Available Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/predict` | POST | Submit a new prediction job |
| `/results/{job_id}` | GET | Get results for a specific job |
| `/status/{job_id}` | GET | Check job status |
| `/health` | GET | API health check |

### Example Results

| Method | Input Structure | Predicted Loop | Complete Structure |
|--------|----------------|----------------|-------------------|
| Simple | [Figure 1] | [Figure 2] | [Figure 3] |
| Modeller | [Figure 4] | [Figure 5] | [Figure 6] |

## Project Structure 📁

```
LoopsAndRibbons/
├── loops_api/
│   ├── utils/
│   │   ├── base_loop_predictor.py    # Base class with common functionality
│   │   ├── simple_predictor.py       # Simple interpolation implementation
│   │   └── modeller_predictor.py     # Modeller-based implementation
│   ├── api/
│   │   ├── routes.py                 # API endpoint definitions
│   │   └── models.py                 # Pydantic models for request/response
│   └── test_loop_predictor.py        # Test suite
├── predictions/                      # Output directory for predictions
├── pyproject.toml                    # Project configuration
└── README.md
```

## Features ✨

### Base Loop Predictor
- Structure loading and parsing (PDB/mmCIF)
- Missing region detection
- Sequence mapping and validation
- Quality score calculation
- Results saving and formatting

### Simple Predictor
- Linear interpolation with random variation
- Fast predictions for quick results
- Basic quality metrics

### Modeller Predictor
- Physics-based loop modeling
- Multiple conformation generation
- Advanced quality assessment
- Complete structure generation

## Quality Metrics 📊

The API calculates quality scores based on:
- Distance from start/end points
- Loop smoothness (angle between segments)
- Bond length consistency

## Limitations ⚠️

### General Limitations
- Maximum input structure size: 100,000 atoms
- Maximum loop length: 30 residues
- Supported file formats: PDB and mmCIF
- No support for multi-chain loop predictions
- No support for non-protein molecules (DNA, RNA, ligands)

### Simple Predictor Limitations
- No consideration of protein physics or energy minimization
- May produce unrealistic conformations for long loops
- Quality scores are based on geometric criteria only
- Not suitable for loops with complex secondary structure

### Modeller Predictor Limitations
- Requires valid Modeller license
- Significantly slower than Simple Predictor
- Memory intensive for large structures
- May fail for very long loops (>20 residues)
- Limited to standard amino acids

## Future Work 🔮

### AlphaFold Integration
- Integration with AlphaFold2 for high-accuracy loop predictions
- Support for multi-chain predictions
- Improved handling of long loops (>30 residues)
- Better prediction of loops with complex secondary structure
- Integration with AlphaFold's confidence metrics

### RareFold Integration
- Implementation of RareFold's specialized loop prediction
- Support for rare loop conformations
- Improved accuracy for challenging loop regions
- Integration with RareFold's quality assessment
- Enhanced handling of non-canonical loop structures

### Planned Improvements
- Hybrid approach combining multiple prediction methods
- Enhanced quality metrics using machine learning
- Support for non-standard amino acids
- Improved handling of post-translational modifications
- Better integration with experimental data
- Real-time prediction progress monitoring
- Batch processing capabilities

## API Documentation 📚

Once the server is running, visit:
- Swagger UI: `http://localhost:8000/docs`
- ReDoc: `http://localhost:8000/redoc`

## License

MIT

## Acknowledgments 🙏

- Modeller team for the excellent modeling software
- Biopython team for the structure handling capabilities
- [Add other acknowledgments]
