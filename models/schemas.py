from pydantic import BaseModel, Field, field_validator, model_validator
from typing import Optional, Literal, List, Any, Dict
from enum import Enum
from datetime import datetime


# ── Enums ─────────────────────────────────────────────────────────────────────


class PipelineSteps(BaseModel):
    drug_likeness: bool = True
    admet: bool = True
    binding_prefilter: bool = True
    docking: bool = True
    retrosynthesis: bool = True


class BindingSiteMode(str, Enum):
    auto = "auto"
    coordinates = "coordinates"
    residues = "residues"


class DockingSpeed(str, Enum):
    fast = "fast"
    balanced = "balanced"
    thorough = "thorough"


class SolubilityFilterMode(str, Enum):
    soluble_only = "soluble_only"
    allow_slightly = "allow_slightly"
    all = "all"


# ── ADMET config models ───────────────────────────────────────────────────────


class ADMETThresholdConfig(BaseModel):
    enabled: bool = True
    cutoff: float
    direction: Literal["above", "below"]
    severity_high: Optional[float] = None
    hard_fail: bool = False
    threshold_str: Optional[str] = None
    implication: Optional[str] = None
    recommendation: Optional[str] = None

    @model_validator(mode="after")
    def validate_threshold_logic(self) -> "ADMETThresholdConfig":
        if self.severity_high is not None:
            if self.direction == "above" and self.severity_high < self.cutoff:
                raise ValueError("For direction='above', severity_high must be >= cutoff")
            if self.direction == "below" and self.severity_high > self.cutoff:
                raise ValueError("For direction='below', severity_high must be <= cutoff")
        return self


class ADMETTuningConfig(BaseModel):
    preset: Literal["balanced", "cns", "oral", "custom"] = "balanced"
    include_full_profile: bool = True

    herg_inhibition: ADMETThresholdConfig = Field(
        default_factory=lambda: ADMETThresholdConfig(
            cutoff=0.50,
            direction="above",
            severity_high=0.85,
            hard_fail=True,
            threshold_str=">0.50",
            implication="cardiac arrhythmia / QT prolongation risk",
            recommendation=(
                "Redesign to remove hERG-liable substructure "
                "(basic amine + aromatic hydrophobic bulk). "
                "Consider bioisostere substitution or pKa reduction."
            ),
        )
    )

    hepatotoxicity: ADMETThresholdConfig = Field(
        default_factory=lambda: ADMETThresholdConfig(
            cutoff=0.50,
            direction="above",
            severity_high=0.85,
            hard_fail=True,
            threshold_str=">0.50",
            implication="liver toxicity risk — potential DILI (drug-induced liver injury)",
            recommendation=(
                "Check for reactive metabolite formation. "
                "Consider metabolic soft-spot analysis; reduce lipophilicity (logP)."
            ),
        )
    )

    caco2_permeability: ADMETThresholdConfig = Field(
        default_factory=lambda: ADMETThresholdConfig(
            cutoff=-5.15,
            direction="below",
            severity_high=-6.5,
            hard_fail=False,
            threshold_str="<-5.15 log cm/s",
            implication="poor intestinal absorption — low oral bioavailability predicted",
            recommendation=(
                "Reduce polar surface area (PSA < 140 Å²) or HBD count. "
                "Consider prodrug strategy or permeation enhancer formulation."
            ),
        )
    )

    oral_bioavailability: ADMETThresholdConfig = Field(
        default_factory=lambda: ADMETThresholdConfig(
            cutoff=0.30,
            direction="below",
            severity_high=0.10,
            hard_fail=False,
            threshold_str="<0.30",
            implication="poor systemic exposure after oral dosing",
            recommendation=(
                "Evaluate first-pass metabolism and P-gp efflux. "
                "Consider IV or alternative delivery route, or structural optimisation."
            ),
        )
    )

    bbb_penetration: ADMETThresholdConfig = Field(
        default_factory=lambda: ADMETThresholdConfig(
            cutoff=0.30,
            direction="below",
            severity_high=None,
            hard_fail=False,
            threshold_str="<0.30",
            implication="limited CNS penetration — may not reach brain target",
            recommendation=(
                "Only relevant if CNS activity is required. "
                "Increase lipophilicity and reduce HBD count to improve CNS exposure."
            ),
        )
    )


# ── Sub-models ────────────────────────────────────────────────────────────────


class BindingSiteCoords(BaseModel):
    x: float
    y: float
    z: float
    box_size: float = 20.0


class BindingSiteResidues(BaseModel):
    chain: str
    residue_start: int
    residue_end: int


class BindingSiteInfo(BaseModel):
    detection_mode: str
    detected_ligand_name: Optional[str] = None
    center_x: float
    center_y: float
    center_z: float
    box_size: float


# ── Job request ───────────────────────────────────────────────────────────────


class JobRequest(BaseModel):
    smiles: str = Field(..., description="SMILES string of the base compound")
    pdb_id: Optional[str] = Field(None, description="RCSB PDB ID (e.g. 1ABC)")
    pdb_content: Optional[str] = Field(None, description="Raw .pdb file content as string")

    num_analogues: int = Field(25, ge=0, le=10000)

    direct_score_only: bool = Field(
        False,
        description="Skip analogue generation; score the input compound directly",
    )

    solubility_filter: SolubilityFilterMode = Field(
        SolubilityFilterMode.all,
        description=(
            "Solubility filter applied at analogue generation time. "
            "'soluble_only' = logS > -3, 'allow_slightly' = logS > -5, 'all' = no filter."
        ),
    )

    toxicity_report_only: bool = Field(
        False,
        description=(
            "Force ADMET toxicity analysis to run and use the ADMET-only shortcut path. "
            "PDB / docking are NEVER required when this is True, even combined with "
            "direct_score_only. Results are always written even if 0 compounds pass ADMET."
        ),
    )

    admet_config: ADMETTuningConfig = Field(
        default_factory=ADMETTuningConfig,
        description=(
            "ADMET decision thresholds and profile-generation settings. "
            "Frontend can override per endpoint."
        ),
    )

    pipeline_steps: PipelineSteps = Field(default_factory=PipelineSteps)
    binding_site_mode: BindingSiteMode = BindingSiteMode.auto
    binding_site_coords: Optional[BindingSiteCoords] = None
    binding_site_residues: Optional[BindingSiteResidues] = None
    docking_speed: DockingSpeed = DockingSpeed.balanced

    max_docking_compounds: int = Field(
        10,
        ge=1,
        le=50,
        description=(
            "Maximum number of compounds ranked by the ML binding pre-filter "
            "to send into AutoDock Vina. Higher values improve hit coverage "
            "but increase docking runtime."
        ),
    )

    mw_min: float = Field(200.0, ge=0, description="Min MW for drug-likeness filter (Da)")
    mw_max: float = Field(500.0, ge=0, description="Max MW for drug-likeness filter (Da)")

    max_lipinski_violations: Optional[int] = Field(
        1,
        description="Max Lipinski violations allowed (1–3). None = ignore filter entirely.",
        ge=1,
        le=3,
    )

    locked_scaffold_smarts: Optional[str] = Field(
        None,
        description=(
            "SMARTS pattern that every generated analogue must contain as a substructure. "
            "None = no scaffold lock. Invalid SMARTS = lock silently disabled with a warning."
        ),
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
        if self.mw_max < self.mw_min:
            raise ValueError("mw_max must be >= mw_min")

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

        if self.toxicity_report_only:
            updated_steps = self.pipeline_steps.model_copy(update={"admet": True})
            object.__setattr__(self, "pipeline_steps", updated_steps)

        return self


# ── Status / steps ────────────────────────────────────────────────────────────


class StepStatus(str, Enum):
    waiting = "waiting"
    running = "running"
    done = "done"
    skipped = "skipped"
    failed = "failed"


class PipelineStepInfo(BaseModel):
    name: str
    status: StepStatus = StepStatus.waiting
    message: Optional[str] = None
    duration_seconds: Optional[float] = None
    progress_current: Optional[int] = None
    progress_total: Optional[int] = None


class JobStatusResponse(BaseModel):
    job_id: str
    status: Literal["queued", "running", "done", "failed"]
    created_at: datetime
    updated_at: datetime
    steps: List[PipelineStepInfo] = Field(default_factory=list)
    error: Optional[str] = None


# ── Pipeline step results ─────────────────────────────────────────────────────


class LipinskiResult(BaseModel):
    passed: bool
    mw: float
    logp: float
    hbd: int
    hba: int
    logs: float
    solubility_class: str


class ADMETFlagDetail(BaseModel):
    property_name: str
    value: float
    threshold: str
    direction: str
    severity: str
    implication: str
    recommendation: str


class ADMETEndpointValue(BaseModel):
    key: str
    label: str
    category: Literal[
        "absorption",
        "distribution",
        "metabolism",
        "excretion",
        "toxicity",
        "physchem",
        "other",
    ] = "other"
    value: float
    display_value: str
    unit: Optional[str] = None
    interpretation: Optional[str] = None
    threshold_applied: Optional[str] = None
    severity: Optional[str] = None
    hard_fail: bool = False
    triggered: bool = False


class ADMETNarrativeBlock(BaseModel):
    title: str
    level: Literal["info", "good", "warning", "danger"] = "info"
    body: str


class ADMETResult(BaseModel):
    passed: bool
    overall_risk: Literal["low", "moderate", "high"] = "low"

    herg_inhibition: float
    caco2_permeability: float
    bbb_penetration: float
    hepatotoxicity: float
    oral_bioavailability: float

    flags: List[ADMETFlagDetail] = Field(default_factory=list)
    flag_summary: List[str] = Field(default_factory=list)
    decision_basis: List[str] = Field(default_factory=list)

    narrative: List[ADMETNarrativeBlock] = Field(default_factory=list)
    endpoint_table: List[ADMETEndpointValue] = Field(default_factory=list)

    extra_properties: Dict[str, float] = Field(default_factory=dict)
    thresholds_used: Dict[str, ADMETThresholdConfig] = Field(default_factory=dict)


class BindingPrefilterResult(BaseModel):
    passed: bool
    predicted_affinity_kcal: float
    confidence: float


class DockingPose(BaseModel):
    rank: int
    affinity_kcal: float
    rmsd_lb: float
    rmsd_ub: float


class DockingResult(BaseModel):
    passed: bool
    best_affinity_kcal: float
    cnn_score: float
    poses: List[DockingPose]


# ── Retrosynthesis ────────────────────────────────────────────────────────────


class RetrosynthesisStep(BaseModel):
    step_number: int
    reaction_smarts: str
    reaction_name: Optional[str] = None
    reaction_type: Optional[str] = None
    starting_materials: List[str]
    product_smiles: Optional[str] = None
    schematic: Optional[str] = None  # "A + B → C" style string
    reagents: List[str] = Field(default_factory=list)
    solvents: List[str] = Field(default_factory=list)
    conditions: Optional[str] = None
    temperature: Optional[str] = None
    duration: Optional[str] = None
    yield_estimate: Optional[str] = None
    protocol_text: Optional[str] = None  # narrative / protocol paragraph
    confidence: float


class RetrosynthesisResult(BaseModel):
    feasible: bool
    num_steps: int
    sa_score: float = 0.0
    route: List[RetrosynthesisStep]
    complexity_score: float
    difficulty_label: str = "unknown"  # "easy" / "moderate" / "hard" / "infeasible"
    synthesis_summary: Optional[str] = None
    estimated_total_yield: Optional[str] = None


class CompoundResult(BaseModel):
    smiles: str
    canonical_smiles: str
    rank: Optional[int] = None
    final_score: Optional[float] = None
    lipinski: Optional[LipinskiResult] = None
    admet: Optional[ADMETResult] = None
    binding_prefilter: Optional[BindingPrefilterResult] = None
    docking: Optional[DockingResult] = None
    retrosynthesis: Optional[RetrosynthesisResult] = None


# ── Job responses ─────────────────────────────────────────────────────────────


class JobResultsResponse(BaseModel):
    job_id: str
    base_smiles: str
    total_analogues_generated: int
    compounds_after_lipinski: int
    compounds_after_admet: int
    compounds_after_prefilter: int
    compounds_docked: int
    final_ranked_compounds: List[CompoundResult]
    binding_site_info: Optional[BindingSiteInfo] = None
    created_at: datetime
    completed_at: datetime


class JobSubmitResponse(BaseModel):
    job_id: str
    message: str
    status_url: str
    results_url: str