import httpx
import asyncio
import os
from pathlib import Path
from typing import Optional


RCSB_PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
JOBS_DIR = Path("storage/jobs")


async def fetch_pdb_by_id(pdb_id: str, job_id: str) -> str:
    """
    Downloads a .pdb file from RCSB by PDB ID.
    Saves it to storage/jobs/{job_id}/receptor_raw.pdb
    Returns the saved file path as a string.
    Raises ValueError if PDB ID not found.
    """
    if not pdb_id or not pdb_id.strip():
        raise ValueError(
            "fetch_pdb_by_id called with empty pdb_id — "
            "this should not happen when toxicity_report_only=True."
        )

    pdb_id = pdb_id.strip().upper()
    url = RCSB_PDB_DOWNLOAD_URL.format(pdb_id=pdb_id)

    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    save_path = job_dir / "receptor_raw.pdb"

    async with httpx.AsyncClient(timeout=30.0) as client:
        response = await client.get(url)

    if response.status_code == 404:
        raise ValueError(
            f"PDB ID '{pdb_id}' not found on RCSB. "
            "Please verify the ID at https://www.rcsb.org"
        )
    if response.status_code != 200:
        raise ValueError(
            f"RCSB returned status {response.status_code} for PDB ID '{pdb_id}'. Try again later."
        )

    content = response.text
    if not content.strip().startswith("HEADER") and "ATOM" not in content:
        raise ValueError(
            f"Downloaded file for '{pdb_id}' does not appear to be a valid PDB file."
        )

    save_path.write_text(content)
    return str(save_path)


async def fetch_pdb_metadata(pdb_id: str) -> dict:
    """
    Fetches basic metadata for a PDB entry (title, organism, resolution)
    using RCSB GraphQL API.
    Returns a dict or empty dict on failure.
    """
    pdb_id = pdb_id.strip().upper()
    query = """
    {
      entry(entry_id: "%s") {
        struct {
          title
        }
        rcsb_entry_info {
          resolution_combined
          polymer_entity_count_protein
        }
        rcsb_entry_container_identifiers {
          entry_id
        }
      }
    }
    """ % pdb_id

    try:
        async with httpx.AsyncClient(timeout=15.0) as client:
            response = await client.post(
                RCSB_GRAPHQL_URL,
                json={"query": query},
                headers={"Content-Type": "application/json"}
            )
        data = response.json()
        entry = data.get("data", {}).get("entry")
        if not entry:
            return {}
        return {
            "pdb_id": pdb_id,
            "title": entry.get("struct", {}).get("title", "Unknown"),
            "resolution_angstrom": (
                entry.get("rcsb_entry_info", {})
                     .get("resolution_combined", [None])[0]
            ),
            "protein_chains": (
                entry.get("rcsb_entry_info", {})
                     .get("polymer_entity_count_protein", "Unknown")
            ),
        }
    except Exception:
        return {}


def save_pdb_from_content(pdb_content: Optional[str], job_id: str) -> str:
    """
    Saves raw .pdb content (from frontend file upload) to disk.
    Returns the saved file path as a string.
    """
    # Guard: pdb_content is None when no file was uploaded and no pdb_id was given.
    # This should never be called in toxicity_report_only mode — but defensive check here.
    if pdb_content is None:
        raise ValueError(
            "save_pdb_from_content called with pdb_content=None — "
            "no PDB file was uploaded. "
            "If you are running toxicity_report_only=True, the orchestrator "
            "should not be calling this function at all. "
            "Check that tox_mode routing in orchestrator.py is correct."
        )

    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    save_path = job_dir / "receptor_raw.pdb"

    content = pdb_content.strip()
    if "ATOM" not in content and "HETATM" not in content:
        raise ValueError(
            "Uploaded file does not appear to be a valid PDB file. "
            "It must contain ATOM or HETATM records."
        )

    save_path.write_text(content)
    return str(save_path)