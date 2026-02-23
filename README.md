# Bio-Engine Sidecar

The Bio-Engine is an asynchronous backend service designed to handle heavy bioinformatics computational tasks for MSAnalyzer. It acts as an integration layer between the main web application and native bioinformatics tools like `tracy`, executing alignments, variant calling, and sequence annotations (HGVS, VEP).

## Features
- **Asynchronous Job Management**: Jobs are executed in isolated background threads.
- **Tracy Wrapper**: Advanced DNA sequence decomposition, basecalling, and alignment using the [tracy](https://github.com/gear-genomics/tracy) C++ framework.
- **Variant Recoder & VEP**: Genomic consequence prediction via the Ensembl REST API with built-in batch handling and fallback strategies.
- **HGVS Notation**: Accurate nomenclature mapping powered by Universal Transcript Archive (UTA).

## System Requirements
- Python >= 3.10
- [Tracy CLI API](https://github.com/gear-genomics/tracy): Must be installed and available in the system `PATH`.
- `samtools` and `bgzip`: Needed for auto-indexing very large reference files (>50Kbp).

## Installation

You can install the dependencies via `pip`:

```bash
pip install -r requirements.txt
```

*(Alternatively, use `pip install .` to install via `pyproject.toml`)*

## Running the Engine

Start the FastAPI server via Uvicorn:

```bash
uvicorn main:app --host 127.0.0.1 --port 8000
```
Or simply:
```bash
python main.py
```

## API Documentation
Once the server is running, the Swagger UI is available at `http://127.0.0.1:8000/docs`, where you can explore and interact with the endpoints.

## Development & Refactoring
This codebase follows PEP8 guidelines enforced by `ruff`. To check style errors, run:
```bash
ruff check .
```
And to automatically format your changes:
```bash
ruff format .
```
