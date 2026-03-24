from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Optional, Literal, List, Any
from enum import Enum
from datetime import datetime


# ── Enums ─────────────────────────────────────────────────────────────────────

class PipelineSteps(BaseModel):
    drug_likeness:    bool = True
    admet:            bool = True
    binding_prefilter: bool = True
    docking:          bool = True
    retrosynthesis:   bool = True


class BindingSiteMode(str, Enum):
    auto        = "auto"
    coordinates = "coordinates"
    residues    = "residues"


class DockingSpeed(str, Enum):
    fast     = "fast"
    balanced = "balanced"
    thorough = "thorough"


class SolubilityFilterMode(str, Enum):
    soluble_only   = "soluble_only"    # logS > -3 (highly soluble + soluble)
    allow_slightly = "allow_slightly"  # logS > -5 (also slightly soluble)
    all            = "all"             # no filter (original behaviour)


# ── Sub-models ────────────────────────────────────────────────────────────────

class BindingSiteCoords(BaseModel):
    x:        float
    y:        float
    z:        float
    box_size: float = 20.0


class BindingSiteResidues(BaseModel):
    chain:         str
    residue_start: int
    residue_end:   int


class BindingSiteInfo(BaseModel):
    detection_mode:       str
    detected_ligand_name: Optional[str] = None
    center_x:             float
    center_y:             float
    center_z:             float
    box_size:             float


# ── Job request ───────────────────────────────────────────────────────────────

class JobRequest(BaseModel):
    smiles:       str           = Field(..., description="SMILES string of the base compound")
    pdb_id:       Optional[str] = Field(None, description="RCSB PDB ID (e.g. 1ABC)")
    pdb_content:  Optional[str] = Field(None, description="Raw .pdb file content as string")

    # ── was Literal[0,10,25,50] — now supports 1–1000 ─────────────────────────
    num_analogues: int = Field(25, ge=0, le=1000)

    direct_score_only: bool = Field(
        False,
        description="Skip analogue generation; score the input compound directly",
    )

    # ── NEW: solubility filter applied during analogue generation ─────────────
    solubility_filter: SolubilityFilterMode = Field(
        SolubilityFilterMode.all,
        description=(
            "Solubility filter applied at analogue generation time. "
            "'soluble_only' = logS > -3, 'allow_slightly' = logS > -5, 'all' = no filter."
        ),
    )

    # ── NEW: skip everything except ADMET — no PDB required ──────────────────
    toxicity_report_only: bool = Field(
        False,
        description=(
            "Run only ADMET toxicity analysis. Skips analogue gen (unless "
            "num_analogues > 0), drug-likeness, docking, and retrosynthesis. "
            "pdb_id / pdb_content are NOT required when this is True."
        ),
    )

    pipeline_steps:        PipelineSteps           = Field(default_factory=PipelineSteps)
    binding_site_mode:     BindingSiteMode          = BindingSiteMode.auto
    binding_site_coords:   Optional[BindingSiteCoords]   = None
    binding_site_residues: Optional[BindingSiteResidues] = None
    docking_speed:         DockingSpeed             = DockingSpeed.balanced

    mw_min: float = Field(200.0, ge=0, description="Min MW for drug-likeness filter (Da)")
    mw_max: float = Field(500.0, ge=0, description="Max MW for drug-likeness filter (Da)")

    max_lipinski_violations: Optional[int] = Field(
        1,
        description="Max Lipinski violations allowed (1–3). None = ignore filter entirely.",
        ge=1,
        le=3,
    )

    @field_validator("smiles")
    @classmethod
    def smiles_not_empty(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("SMILES string cannot be empty")
        return v

    @field_validator("pdb_id")
    @classmethod
    def pdb_id_format(cls, v: Optional[str]) -> Optional[str]:
        if v is not None:
            v = v.strip().upper()
            if len(v) != 4 or not v.isalnum():
                raise ValueError("PDB ID must be exactly 4 alphanumeric characters")
        return v

    @model_validator(mode="after")
    def cross_field_checks(self) -> "JobRequest":
        # PDB only required when docking is actually going to run
        needs_pdb = (
            not self.toxicity_report_only
            and self.pipeline_steps.docking
        )
        if needs_pdb and self.pdb_id is None and self.pdb_content is None:
            raise ValueError(
                "Either pdb_id or pdb_content must be provided "
                "(not required when toxicity_report_only=True or docking is disabled)"
            )
        if self.binding_site_mode == BindingSiteMode.coordinates and self.binding_site_coords is None:
            raise ValueError("binding_site_coords required when mode is 'coordinates'")
        if self.binding_site_mode == BindingSiteMode.residues and self.binding_site_residues is None:
            raise ValueError("binding_site_residues required when mode is 'residues'")
        if self.direct_score_only:
            object.__setattr__(self, "num_analogues", 0)
        return self


# ── Status / steps ────────────────────────────────────────────────────────────

class StepStatus(str, Enum):
    waiting = "waiting"
    running = "running"
    done    = "done"
    skipped = "skipped"
    failed  = "failed"


class PipelineStepInfo(BaseModel):
    name:             str
    status:           StepStatus    = StepStatus.waiting
    message:          Optional[str] = None
    duration_seconds: Optional[float] = None


class JobStatusResponse(BaseModel):
    job_id:     str
    status:     Literal["queued", "running", "done", "failed"]
    created_at: datetime
    updated_at: datetime
    steps:      List[PipelineStepInfo] = Field(default_factory=list)
    error:      Optional[str]          = None


# ── Pipeline step results ─────────────────────────────────────────────────────

class LipinskiResult(BaseModel):
    passed:          bool
    mw:              float
    logp:            float
    hbd:             int
    hba:             int
    logs:            float
    solubility_class: str


# ── NEW: structured ADMET flag with full context ──────────────────────────────

class ADMETFlagDetail(BaseModel):
    property_name:  str    # e.g. "hERG Inhibition"
    value:          float  # e.g. 0.87
    threshold:      str    # e.g. ">0.50"
    direction:      str    # "above" | "below"
    severity:       str    # "high" | "moderate" | "low"
    implication:    str    # e.g. "cardiac arrhythmia / QT prolongation risk"
    recommendation: str    # redesign suggestion


class ADMETResult(BaseModel):
    passed:               bool
    herg_inhibition:      float
    caco2_permeability:   float
    bbb_penetration:      float
    hepatotoxicity:       float
    oral_bioavailability: float
    # ← was List[str] — now structured objects for the report
    flags:        List[ADMETFlagDetail] = []
    # ← plain-English one-liners for the UI badge list
    flag_summary: List[str]            = []


class BindingPrefilterResult(BaseModel):
    passed:                  bool
    predicted_affinity_kcal: float
    confidence:              float


class DockingPose(BaseModel):
    rank:        int
    affinity_kcal: float
    rmsd_lb:     float
    rmsd_ub:     float


class DockingResult(BaseModel):
    passed:            bool
    best_affinity_kcal: float
    cnn_score:         float
    poses:             List[DockingPose]


class RetrosynthesisStep(BaseModel):
    step_number:        int
    reaction_smarts:    str
    starting_materials: List[str]
    confidence:         float


class RetrosynthesisResult(BaseModel):
    feasible:         bool
    num_steps:        int
    route:            List[RetrosynthesisStep]
    complexity_score: float


class CompoundResult(BaseModel):
    smiles:           str
    canonical_smiles: str
    rank:             Optional[int]   = None
    final_score:      Optional[float] = None
    lipinski:         Optional[LipinskiResult]         = None
    admet:            Optional[ADMETResult]            = None
    binding_prefilter: Optional[BindingPrefilterResult] = None
    docking:          Optional[DockingResult]          = None
    retrosynthesis:   Optional[RetrosynthesisResult]   = None


# ── Job responses ─────────────────────────────────────────────────────────────

class JobResultsResponse(BaseModel):
    job_id:                    str
    base_smiles:               str
    total_analogues_generated: int
    compounds_after_lipinski:  int
    compounds_after_admet:     int
    compounds_after_prefilter: int
    compounds_docked:          int
    final_ranked_compounds:    List[CompoundResult]
    binding_site_info:         Optional[BindingSiteInfo] = None
    created_at:                datetime
    completed_at:              datetime


class JobSubmitResponse(BaseModel):
    job_id:      str
    message:     str
    status_url:  str
    results_url: str
