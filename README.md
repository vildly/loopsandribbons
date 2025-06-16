# LoopsAndRibbons API 🧬

A RESTful API service for predicting missing loops in protein structures using both simple interpolation and Modeller-based approaches.

## Overview 🔬

This project provides tools for identifying and predicting missing regions (loops) in protein structures. It offers two prediction methods:
1. **Simple Predictor**: Uses linear interpolation with random variation for quick predictions
2. **Modeller Predictor**: Uses Modeller for more accurate, physics-based predictions

---

## Current Implementation ✅

- **Core Utilities and Classes** (in `loops_api/utils/`):
  - `base_loop_predictor.py`: Base class for predictors, structure parsing, sequence mapping, and region detection
  - `simple_predictor.py`: Simple interpolation-based loop predictor
  - `modeller_predictor.py`: Modeller-based loop predictor
  - `prediction_result.py`: Result writing and summary utilities
  - `loop_assembler.py`: Assembles final structures, handles renumbering, and generates Ramachandran plots
- **Test Script**
  - `test_loop_predictor.py`: End-to-end tests for predictors and structure assembly
- **Results and Examples**
  - Example output images and PDBs in `predictions/`
  - Example Ramachandran plots and structure visualizations
- **README**
  - Example results table with images
  - Project structure and feature overview

---

## Not Yet Implemented / Future Work 🚧

- **API Endpoints**
  - The FastAPI endpoints described below are **not yet implemented**. They are planned for a future release.
  - No running REST API server is provided yet
- **Dockerization**
  - No Dockerfile or containerization is provided yet
- **AlphaFold & RareFold Integration**
  - No AlphaFold or RareFold support yet (see Future Work below)
- **Production Deployment**
  - No deployment scripts or cloud setup
- **Batch Processing & Monitoring**
  - No real-time progress or batch job support

---

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

## API Endpoints 🔌 *(Planned)*

> **Note:** The following endpoints are **not yet implemented**. This project does not currently provide any API or server. Making this a real API is planned for future work.

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

### Planned Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/predict` | POST | Submit a new prediction job |
| `/results/{job_id}` | GET | Get results for a specific job |
| `/status/{job_id}` | GET | Check job status |
| `/health` | GET | API health check |

---

# Example Results

| Missing Region Models | Structure & Missing | Ramachandran Plot |
|----------------------|---------------------|-------------------|
| ![](predictions/3idp_missing_region_models.png) | ![](predictions/2idp_and_missing.png) | ![](predictions/3IDP_assembled_ramachandran.png) |

<!-- Add more examples as needed -->

---
## Example Prediction Summary

A `prediction_summary.json` file is generated in each prediction directory. Example:

```json
{
  "timestamp": "20250616_120450",
  "region_metadata": {
    "chain_id": "B",
    "start_res": 597,
    "end_res": 614,
    "length": 16,
    "missing_sequence": "ATEKSRWSGSHQFEQL",
    "full_chain_sequence": "HHHHHHDRNRMKTLGRRDSSDDWEIPDGQITVGQRIGSGSFGTVYKGKWHGDVAVKMLNVTAPTPQQLQAFKNEVGVLRKTRHVNILLFMGYSTKPQLAIVTQWCEGSSLYHHLHIIETKFEMIKLIDIARQTAQGMDYLHAKSIIHRDLKSNNIFLHEDLTVKIGDFGLATEKSRWSGSHQFEQLSGSILWMAPEVIRMQDKNPYSFQSDVYAFGIVLYELMTGQLPYSNINNRDQIIFMVGRGYLSPDLSKVRSNCPKAMKRLMAECLKKKRDERPLFPQILASIELLARSLPKIHRS"
  },
  "num_models": 1,
  "models": [
    {
      "model_number": 1,
      "quality_score": -31798.841796875,
      "ga341_score": [
        1.0,
        0.30858367681503296,
        -340.16363525390625,
        -13.651281356811523,
        -9.861736297607422,
        -8.024914741516113,
        -6.65782356262207,
        -11.363975524902344
      ],
      "model_file": "model.B99990001.pdb"
    }
  ]
}
```

## Project Structure 📁

```
LoopsAndRibbons/
├── loops_api/
│   ├── utils/
│   │   ├── base_loop_predictor.py    # Base class with common functionality
│   │   ├── simple_predictor.py       # Simple interpolation implementation
│   │   └── modeller_predictor.py     # Modeller-based implementation
│   │   ├── prediction_result.py      # Result writing utilities
│   │   └── loop_assembler.py         # Structure assembly and analysis
│   ├── api/
│   │   ├── routes.py                 # (Placeholder) API endpoint definitions
│   │   └── models.py                 # (Placeholder) Pydantic models for request/response
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

### Loop Assembler
- Robust structure merging and renumbering
- Ramachandran plot generation
- Output in PDB format

## Quality Metrics 📊

The quality of predicted loops and models is assessed using:

**For Modeller-based predictions:**
- **DOPE Score**: Discrete Optimized Protein Energy, a statistical potential used by Modeller to evaluate model quality
- **GA341 Score**: (if available) A Modeller score for model reliability
- **Other Model-Specific Scores**: Any additional scores provided by the prediction method (e.g., energy, confidence)
- All scores are included in the `prediction_summary.json` for each prediction

**For Simple Predictor:**
- **Distance from start/end points**
- **Loop smoothness** (angle between segments)
- **Bond length consistency**

For simple predictors, geometric criteria are primary. For Modeller-based predictions, the above model-specific scores are primary.

## Limitations ⚠️

### General Limitations
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

### AlphaFold Integration *(Planned)*
- Integration with AlphaFold2 for high-accuracy loop predictions
- Support for multi-chain predictions
- Improved handling of long loops (>30 residues)
- Better prediction of loops with complex secondary structure
- Integration with AlphaFold's confidence metrics

### RareFold Integration *(Planned)*
- Implementation of RareFold's specialized loop prediction
- Support for rare loop conformations
- Improved accuracy for challenging loop regions
- Integration with RareFold's quality assessment
- Enhanced handling of non-canonical loop structures

### Other Planned Improvements
- **REST API endpoints and server** (full implementation)
- **Dockerization** (containerized deployment)
- **CI/CD pipeline** (automated testing and deployment)
- **Hosting/Deployment** (host the service on a public or private server/cloud)
- Hybrid approach combining multiple prediction methods
- Enhanced quality metrics using machine learning
- Support for non-standard amino acids
- Improved handling of post-translational modifications
- Better integration with experimental data
- Real-time prediction progress monitoring
- Batch processing capabilities

## API Documentation 📚 *(Planned)*

> **Note:** There is currently no API documentation because there is no API yet. Making this project into a real API with documentation is future work.

Once the server is running, visit:
- Swagger UI: `http://localhost:8000/docs`
- ReDoc: `http://localhost:8000/redoc`




## License

MIT