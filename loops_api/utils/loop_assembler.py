from typing import Dict, Tuple, List, Optional
from pathlib import Path
import logging
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBParser, MMCIFParser, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
import numpy as np

from .base_loop_predictor import LoopPredictor

from ramachandraw.parser import get_phi_psi
from ramachandraw.utils import plot

class LoopAssembler:
    """Handles assembly of loop predictions into complete structures with correct numbering"""
    
    def __init__(self, predictor: LoopPredictor):
        """Initialize with the predictor that has the sequence mappings
        
        Args:
            predictor: A LoopPredictor instance that has built sequence mappings
        """
        self.predictor = predictor
        self.original_structure = predictor.structure
        self.logger = logging.getLogger("loop_assembler")
        
    def assemble_structure(self, model_structure: Structure, region_chain_id: str, 
                          first_auth_seq_id: int, last_auth_seq_id: int, 
                          ramachandran_png_path: Optional[str] = None) -> Structure:
        """Assemble a new structure combining original and modeled parts, and optionally save a Ramachandran plot.
        
        Args:
            model_structure: Structure containing the modeled chain
            region_chain_id: Chain ID of the region being modeled
            first_auth_seq_id: First auth_seq_id in the region
            last_auth_seq_id: Last auth_seq_id in the region
            ramachandran_png_path: If provided, save a Ramachandran plot PNG to this path
        Returns:
            A new Structure object containing the assembled structure
        """
        print(f"[LoopAssembler] Assembling structure for chain {region_chain_id}, region {first_auth_seq_id}-{last_auth_seq_id}")
        # Create new structure
        new_structure = Structure.Structure("assembled")
        new_model = Model.Model(0)
        new_structure.add(new_model)

        # Copy all chains from the original structure
        original_chains = {}
        for model in self.original_structure:
            for chain in model:
                print(f"[LoopAssembler] Copying chain {chain.id} from original structure with {len(list(chain))} residues")
                # Deep copy all chains
                new_chain = Chain.Chain(chain.id)
                for residue in chain:
                    new_chain.add(residue.copy())
                new_model.add(new_chain)
                original_chains[chain.id] = new_chain

        # Identify the modeled chain in the MODELLER output (usually only one chain, or chain A)
        modeled_chain = None
        for model in model_structure:
            chains = list(model)
            print(f"[LoopAssembler] MODELLER output chains: {[c.id for c in chains]}")
            if len(chains) == 1:
                modeled_chain = chains[0]
                print(f"[LoopAssembler] Using only chain in model: {modeled_chain.id}")
                break
            for chain in chains:
                if chain.id in ("A", "", region_chain_id):
                    modeled_chain = chain
                    print(f"[LoopAssembler] Selected modeled chain: {chain.id}")
                    break
            if modeled_chain:
                break
        if not modeled_chain:
            raise ValueError(f"Modeled chain not found in model structure (checked for chain A, blank, or {region_chain_id})")

        # Patch only the missing region residues in the modeled chain into the correct chain in the new structure
        chain_to_patch = original_chains.get(region_chain_id)
        if not chain_to_patch:
            raise ValueError(f"Chain {region_chain_id} not found in original structure")
        # Remove residues in the missing region
        residues_to_remove = [res for res in chain_to_patch if is_aa(res) and first_auth_seq_id <= res.id[1] <= last_auth_seq_id]
        print(f"[LoopAssembler] Removing {len(residues_to_remove)} residues from chain {region_chain_id} in region {first_auth_seq_id}-{last_auth_seq_id}")
        for res in residues_to_remove:
            print(f"[LoopAssembler] Removing residue: {res.resname} {res.id}")
            chain_to_patch.detach_child(res.id)

        # --- NEW LOGIC: Map modeller_rel_res_id to label_seq_id to auth_seq_id/auth_icode ---
        # Find the first label_seq_id for the missing region in the mapping
        # This is the minimum label_seq_id that maps to an auth_seq_id in the missing region
        label_seq_ids_in_region = [label_seq_id for (chain, label_seq_id), (auth_seq_id, _) in self.predictor.label_seq_id_to_auth_id_and_icode.items()
                                   if chain == region_chain_id and first_auth_seq_id <= auth_seq_id <= last_auth_seq_id]
        if not label_seq_ids_in_region:
            raise ValueError(f"No label_seq_ids found for chain {region_chain_id} in region {first_auth_seq_id}-{last_auth_seq_id}")
        first_resolved_label_id_of_segment = min(label_seq_ids_in_region)
        print(f"[LoopAssembler] first_resolved_label_id_of_segment: {first_resolved_label_id_of_segment}")

        # Add/replace residues from the modeled chain, renumbering as needed
        for res in modeled_chain:
            if is_aa(res):
                modeller_rel_res_id = res.get_id()[1]  # 1-based in MODELLER output
                corresponding_label_seq_id = first_resolved_label_id_of_segment + (modeller_rel_res_id - 1)
                try:
                    auth_seq_id, auth_icode = self.predictor.get_auth_seq_id(region_chain_id, corresponding_label_seq_id)
                    if first_auth_seq_id <= auth_seq_id <= last_auth_seq_id:
                        print(f"[LoopAssembler] Adding residue: modeller_rel_res_id={modeller_rel_res_id}, label_seq_id={corresponding_label_seq_id}, auth_seq_id={auth_seq_id}, auth_icode={auth_icode}, resname={res.resname}")
                        # Create new residue with correct numbering
                        new_residue = Residue.Residue((" ", auth_seq_id, auth_icode), res.resname, res.segid)
                        for atom in res:
                            new_residue.add(atom.copy())
                        chain_to_patch.add(new_residue)
                except ValueError as e:
                    print(f"[LoopAssembler] Skipping residue {modeller_rel_res_id} (label_seq_id {corresponding_label_seq_id}): {e}")

        print(f"[LoopAssembler] Finished assembling structure.")

        # Optionally generate and save Ramachandran plot
        if ramachandran_png_path is not None:
            import tempfile
            io = PDBIO()
            with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
                io.set_structure(new_structure)
                io.save(tmp.name)
                tmp_pdb_path = tmp.name
            print(f"[LoopAssembler] Generating Ramachandran plot for {tmp_pdb_path}")
            try:
                plot(tmp_pdb_path, outpath=ramachandran_png_path)
                print(f"[LoopAssembler] Ramachandran plot saved to {ramachandran_png_path}")
            except Exception as e:
                print(f"[LoopAssembler] Error generating Ramachandran plot: {e}")

        return new_structure 

    def save_ramachandran_plot(self, structure: Structure, output_png: str):
        """Generate and save a Ramachandran plot for the given structure using ramachandraw."""
        # Save structure to a temporary PDB file
        import tempfile
        io = PDBIO()
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            io.set_structure(structure)
            io.save(tmp.name)
            tmp_pdb_path = tmp.name
        print(f"[LoopAssembler] Generating Ramachandran plot for {tmp_pdb_path}")
        # Generate and save the plot
        try:
            plot(tmp_pdb_path, cmap="viridis", alpha=0.75, dpi=100, save=True, show=False, filename=output_png)
            print(f"[LoopAssembler] Ramachandran plot saved to {output_png}")
        except Exception as e:
            print(f"[LoopAssembler] Error generating Ramachandran plot: {e}") 