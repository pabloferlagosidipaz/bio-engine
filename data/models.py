from __future__ import annotations

import uuid
from enum import Enum
from typing import Any

from pydantic import BaseModel, ConfigDict, Field, model_validator


class Comment(BaseModel):
    id: str
    text: str
    author: str
    created_at: str


class TracyConfig(BaseModel):
    model_config = ConfigDict(extra='ignore')

    # Generic options
    pratio: float = 0.33
    kmer: int = 15
    support: int = 3
    maxindel: int = 1000
    callVariants: bool = True
    vep: bool = False

    # Alignment options
    gapopen: int = -10
    gapext: int = -4
    match: int = 3
    mismatch: int = -5

    # Trimming options
    trim: int = 0
    trimLeft: int = 50
    trimRight: int = 50

    # Metadata & Formatting
    annotate: str | None = None
    linelimit: int = 10000

class HGVSConfig(BaseModel):
    transcript: str | None = None
    gene: str | None = None
    assembly: str = "GRCh38"
    auto_hgvs: bool = False
    auto_vep: bool = False

class JobStatus(str, Enum):
    CREATED = "created"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

class JobReference(BaseModel):
    type: str  # "ncbi" or "file"
    value: str # accession code or file path

class JobRead(BaseModel):
    id: str = Field(default_factory=lambda: str(uuid.uuid4())[:8])
    file: str
    trimLeft: int | None = None
    trimRight: int | None = None

class JobPatient(BaseModel):
    id: str
    name: str
    reads: list[JobRead | dict[str, Any]] # Backward compatibility for raw dicts

    @model_validator(mode='before')
    @classmethod
    def parse_reads(cls, data: Any) -> Any:
        # If reads are strings (legacy), convert to dicts
        if isinstance(data, dict) and 'reads' in data:
            new_reads = []
            for r in data['reads']:
                if isinstance(r, str):
                    new_reads.append({'file': r, 'id': str(uuid.uuid4())[:8]})
                else:
                    new_reads.append(r)
            data['reads'] = new_reads
        return data

class Job(BaseModel):
    id: str
    name: str
    status: JobStatus = JobStatus.CREATED
    progress: int = 0
    status_message: str | None = None
    created_at: str
    updated_at: str
    reference: dict[str, str] # e.g. {"type": "ncbi", "value": "NM_000000"}
    patients: list[dict[str, Any]] # List of patient dicts
    results: list[dict[str, Any]] | None = None # List of alignment results
    reference_sequence: str | None = None
    features: list[dict[str, Any]] | None = None # List of gene features (exons, CDS, etc.)
    app_version: str | None = None
    config: TracyConfig | None = None
    hgvs_config: HGVSConfig | None = None

    # Use Field(default_factory=...) to avoid mutable default issues
    comments: dict[str, list[Comment]] = Field(default_factory=dict)
    hgvs_alternatives: dict[str, list[str]] = Field(default_factory=dict)
    vep_annotations: dict[str, VariantAnnotation] = Field(default_factory=dict)

class CreateJobRequest(BaseModel):
    name: str
    reference: dict[str, str]
    patients: list[dict[str, Any]]
    app_version: str | None = None
    config: TracyConfig | None = None
    hgvs_config: HGVSConfig | None = None

class RenameJobRequest(BaseModel):
    name: str

class UpdateJobRequest(BaseModel):
    name: str
    reference: dict[str, str]
    patients: list[dict[str, Any]]
    config: TracyConfig | None = None
    hgvs_config: HGVSConfig | None = None

class AddCommentRequest(BaseModel):
    variant_key: str
    text: str
    author: str

class AddHGVSAlternativesRequest(BaseModel):
    principal_hgvs: str
    alternatives: list[str]

class FilePaths(BaseModel):
    ref_path: str
    sample_path: str
    config: TracyConfig | None = None
    hgvs_config: HGVSConfig | None = None

class VariantData(BaseModel):
    columns: list[str]
    rows: list[list[Any]]

class VariantAnnotation(BaseModel):
    model_config = ConfigDict(extra='allow')

    gene_symbol: str | None = None
    impact: str | None = None
    consequence: str | None = None
    hgvs_c: str | None = None
    hgvs_p: str | None = None
    sift: str | None = None
    polyphen: str | None = None
    # Enrich as needed from VEP

class AlignmentMetadata(BaseModel):
    refStart: int
    refForward: int

class ConsensusAlignItem(BaseModel):
    refPos: int
    sangerPos: int | None = None
    sangerPos1: list[int] = []
    sangerPos2: list[int] = []
    alt1: list[str]
    alt2: list[str]
    cons: list[str]

class TraceData(BaseModel):
    traceA: list[int]
    traceC: list[int]
    traceG: list[int]
    traceT: list[int]
    peakLocations: list[int]

class AlignmentResponse(BaseModel):
    model_config = ConfigDict(extra='allow')

    variants: VariantData
    consensusAlign: dict[str, ConsensusAlignItem] | None = None
    readSeqConsensus: list[str] | None = None
    readSeqConsensusComplementary: list[str] | None = None
    readSeqRef: list[str] | None = None
    alignment: AlignmentMetadata
    trace: TraceData
    baseCount: int | None = None

    @model_validator(mode='before')
    @classmethod
    def map_tracy_fields(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data

        # Pull standard Tracy fields into structured sub-models if missing
        if 'alignment' not in data and 'ref1pos' in data:
            is_reverse = data.get("ms_analyzer", {}).get("original_orientation", "forward") == "reverse"
            data['alignment'] = {
                "refStart": int(data.get("ref1pos", 0)) - 1,
                "refForward": int(not is_reverse)
            }

        if 'trace' not in data and 'peakA' in data:
            data['trace'] = {
                "traceA": data.get("peakA", []),
                "traceC": data.get("peakC", []),
                "traceG": data.get("peakG", []),
                "traceT": data.get("peakT", []),
                "peakLocations": data.get("basecallPos", [])
            }

        # Fix consensusAlign items missing refPos (robustness)
        if 'consensusAlign' in data and isinstance(data['consensusAlign'], dict):
            for ref_p, item in data['consensusAlign'].items():
                if isinstance(item, dict) and 'refPos' not in item:
                    try:
                        item['refPos'] = int(ref_p)
                    except (ValueError, TypeError):
                        pass

        if 'readSeqRef' in data and isinstance(data['readSeqRef'], str):
            data['readSeqRef'] = list(data['readSeqRef'])

        return data

class HGVSRequest(BaseModel):
    transcript: str
    assembly: str = "GRCh38"
    pos: int
    ref: str
    alt: str
