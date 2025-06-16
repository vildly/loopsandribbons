from typing import Dict, Tuple, List, Optional
from pathlib import Path
import logging
from Bio.PDB import Structure, Model, Chain, Residue, Atom, PDBParser, MMCIFParser, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
import numpy as np

class LoopAssembler:
    """Handles assembly of loop predictions into complete structures with correct numbering"""
    
    def __init__(self, original_structure: Structure):
        """Initialize with the original structure"""
        self.original_structure = original_structure
        self.label_seq_id_to_auth_id_and_icode: Dict[Tuple[str, int], Tuple[int, str]] = {}
        self.auth_id_to_label_seq_id: Dict[Tuple[str, int], int] = {}
        self.logger = logging.getLogger("loop_assembler")
        
    def build_sequence_mappings(self, chain_id: str, first_auth_seq_id: int, last_auth_seq_id: int) -> None:
        """Build mappings between label_seq_id and auth_seq_id for a chain
        
        Args:
            chain_id: The chain ID to build mappings for
            first_auth_seq_id: First auth_seq_id in the chain
            last_auth_seq_id: Last auth_seq_id in the chain
        """
        # Get the chain from the original structure
        chain = None
        for model in self.original_structure:
            if chain_id in model:
                chain = model[chain_id]
                break
                
        if not chain:
            raise ValueError(f"Chain {chain_id} not found in original structure")
            
        # Build the mapping for observed residues
        current_label_seq_id = 1
        for residue in sorted(chain.get_residues(), key=lambda r: r.id[1]):
            if is_aa(residue):
                auth_seq_id = residue.id[1]
                auth_icode = residue.id[2]
                
                # Store both mappings
                self.label_seq_id_to_auth_id_and_icode[(chain_id, current_label_seq_id)] = (auth_seq_id, auth_icode)
                self.auth_id_to_label_seq_id[(chain_id, auth_seq_id)] = current_label_seq_id
                
                current_label_seq_id += 1
                
        # For missing residues within the range, extend the mapping linearly
        for auth_seq_id in range(first_auth_seq_id, last_auth_seq_id + 1):
            if (chain_id, auth_seq_id) not in self.auth_id_to_label_seq_id:
                # Find the next available label_seq_id
                while (chain_id, current_label_seq_id) in self.label_seq_id_to_auth_id_and_icode:
                    current_label_seq_id += 1
                    
                # Add the mapping for this missing residue
                self.label_seq_id_to_auth_id_and_icode[(chain_id, current_label_seq_id)] = (auth_seq_id, ' ')
                self.auth_id_to_label_seq_id[(chain_id, auth_seq_id)] = current_label_seq_id
                current_label_seq_id += 1
                
        self.logger.info(f"Built sequence mappings for chain {chain_id}")
        self.logger.debug(f"Number of mappings: {len(self.label_seq_id_to_auth_id_and_icode)}")
        
    def assemble_structure(self, model_structure: Structure, region_chain_id: str, 
                          first_auth_seq_id: int, last_auth_seq_id: int) -> Structure:
        """Assemble a new structure combining original and modeled parts
        
        Args:
            model_structure: Structure containing the modeled chain
            region_chain_id: Chain ID of the region being modeled
            first_auth_seq_id: First auth_seq_id in the region
            last_auth_seq_id: Last auth_seq_id in the region
            
        Returns:
            A new Structure object containing the assembled structure
        """
        # Create new structure
        new_structure = Structure.Structure("assembled")
        new_model = Model.Model(0)
        new_structure.add(new_model)
        
        # Process each chain in the original structure
        for model in self.original_structure:
            for chain_id, chain in model.child_dict.items():
                if chain_id != region_chain_id:
                    # For non-modeled chains, copy everything
                    new_chain = Chain.Chain(chain_id)
                    for residue in chain:
                        new_chain.add(residue.copy())
                    new_model.add(new_chain)
                else:
                    # For the modeled chain, only copy HETATMs and residues outside the region
                    new_chain = Chain.Chain(chain_id)
                    for residue in chain:
                        if not is_aa(residue) or residue.id[1] < first_auth_seq_id or residue.id[1] > last_auth_seq_id:
                            new_chain.add(residue.copy())
                    new_model.add(new_chain)
        
        # Get the modeled chain and renumber it
        modeled_chain = None
        for model in model_structure:
            if region_chain_id in model:
                modeled_chain = model[region_chain_id]
                break
                
        if not modeled_chain:
            raise ValueError(f"Modeled chain {region_chain_id} not found in model structure")
            
        # Create a new chain for the modeled region
        new_modeled_chain = Chain.Chain(region_chain_id)
        
        # Add and renumber the modeled residues
        for residue in modeled_chain:
            if is_aa(residue):
                # Get the label_seq_id for this position
                label_seq_id = residue.id[1]  # Modeller uses 1-based numbering
                
                # Get the corresponding auth_seq_id and icode
                if (region_chain_id, label_seq_id) in self.label_seq_id_to_auth_id_and_icode:
                    auth_seq_id, auth_icode = self.label_seq_id_to_auth_id_and_icode[(region_chain_id, label_seq_id)]
                    
                    # Create new residue with correct numbering
                    new_residue = Residue.Residue((auth_seq_id, auth_icode, ' '), residue.resname, residue.segid)
                    
                    # Copy all atoms
                    for atom in residue:
                        new_residue.add(atom.copy())
                        
                    new_modeled_chain.add(new_residue)
                    
        # Add the renumbered modeled chain
        new_model.add(new_modeled_chain)
        
        return new_structure 