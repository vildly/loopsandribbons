import numpy as np
from typing import List, Dict, Tuple, Optional
from pathlib import Path
import json
from Bio.PDB import PDBParser, MMCIFParser, Structure, Model, Chain, Residue, Atom, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import logging
from .missing_region import MissingRegion
from .prediction_result import ResultWriter, FileResultWriter

class LoopPredictor:
    """Base class for loop prediction methods"""
    
    def __init__(self, pdb_file_path: str, result_writer: Optional[ResultWriter] = None):
        self.pdb_file_path = pdb_file_path
        self.structure = None
        self.cif_dict = {}
        self.full_chain_sequences: Dict[str, str] = {}
        self.observed_residues_by_label_seq_id: Dict[Tuple[str, int], Tuple[int, str, Residue.Residue]] = {}
        self.label_seq_id_to_full_seq_idx: Dict[Tuple[str, int], int] = {}
        self.original_flanking_residues: Dict[Tuple[str, int, str], Residue.Residue] = {}
        self.auth_id_to_label_seq_id: Dict[Tuple[str, int], int] = {}
        self.result_writer = result_writer or FileResultWriter()

    def _get_chain_sequence_from_mmcif(self, chain_id: str) -> Optional[str]:
        """Get the sequence from mmCIF _entity_poly_seq or _struct_ref_seq records"""
        if not isinstance(self.structure.header, dict):
            return None
            
        # Try _entity_poly_seq first (most reliable)
        if '_entity_poly_seq' in self.structure.header:
            for record in self.structure.header['_entity_poly_seq']:
                if record.get('pdbx_strand_id') == chain_id:
                    try:
                        return ''.join(seq1(res) for res in record['mon_id'])
                    except:
                        continue
        
        # Fallback to _struct_ref_seq
        if '_struct_ref_seq' in self.structure.header:
            for record in self.structure.header['_struct_ref_seq']:
                if record.get('pdbx_strand_id') == chain_id:
                    try:
                        return ''.join(seq1(res) for res in record['mon_id'])
                    except:
                        continue
        
        return None

    def _create_residue_mapping(self, chain) -> Dict[int, int]:
        """Create a mapping between PDB residue numbers and sequence indices"""
        mapping = {}
        seq_idx = 0
        
        # Get the full sequence for this chain
        chain_id = chain.id
        if chain_id not in self.full_chain_sequences:
            return mapping
            
        # Create mapping by iterating through residues in order
        for residue in chain:
            if is_aa(residue):
                mapping[residue.id[1]] = seq_idx
                seq_idx += 1
                
        return mapping

    def _extract_chain_sequences_from_mmcif(self):
        """Extract full chain sequences for each chain from mmCIF _entity_poly and _entity_poly_seq"""
        if not hasattr(self.structure, 'header') or not isinstance(self.structure.header, dict):
            return
        header = self.structure.header
        # Map entity_id to chain_ids
        entity_to_chains = {}
        if '_entity_poly' in header:
            for record in header['_entity_poly']:
                entity_id = record['entity_id'] if 'entity_id' in record else record.get('_row', None)
                if not entity_id:
                    continue
                chain_ids = record['pdbx_strand_id'].replace(' ', '').split(',')
                entity_to_chains[entity_id] = chain_ids
                # Get the one-letter sequence for this entity
                seq = record['pdbx_seq_one_letter_code'].replace('\n', '').replace('\r', '').replace(' ', '').replace(';', '')
                for chain_id in chain_ids:
                    self.full_chain_sequences[chain_id] = seq
        # Store mapping from entity_id to (num, mon_id) for residue mapping
        self.entity_poly_seq = {}
        if '_entity_poly_seq' in header:
            for record in header['_entity_poly_seq']:
                entity_id = str(record['entity_id'])
                num = int(record['num'])
                mon_id = record['mon_id']
                if entity_id not in self.entity_poly_seq:
                    self.entity_poly_seq[entity_id] = []
                self.entity_poly_seq[entity_id].append((num, mon_id))
        self.chain_to_entity = {}
        for entity_id, chain_ids in entity_to_chains.items():
            for chain_id in chain_ids:
                self.chain_to_entity[chain_id] = entity_id

    def _get_chain_sequence(self, chain) -> str:
        chain_id = chain.id
        # If already cached, return
        if chain_id in self.full_chain_sequences:
            return self.full_chain_sequences[chain_id]
        # For mmCIF, extract from header
        if self.pdb_file_path.endswith('.cif'):
            self._extract_chain_sequences_from_mmcif()
            return self.full_chain_sequences.get(chain_id, '')
        # Fallback to ATOM/HETATM records
        sequence = []
        for residue in chain:
            try:
                # Try to convert to one-letter code, fallback to X if not possible
                sequence.append(seq1(residue.resname))
            except:
                sequence.append('X')
        sequence = ''.join(sequence)
        self.full_chain_sequences[chain_id] = sequence
        return sequence

    def _get_missing_sequence(self, chain_id: str, start_res: int, end_res: int) -> str:
        """Get the sequence of missing residues using residue mapping"""
        if chain_id not in self.full_chain_sequences:
            # Fallback to X's if we can't get the actual sequence
            return 'X' * (end_res - start_res - 1)
            
        full_sequence = self.full_chain_sequences[chain_id]
        
        # Convert PDB residue numbers to label sequence IDs
        start_label_seq_id = None
        end_label_seq_id = None
        
        # Look up the label_seq_ids using auth_seq_ids
        for (cid, auth_id), label_seq_id in self.auth_id_to_label_seq_id.items():
            if cid == chain_id:
                if auth_id == start_res:
                    start_label_seq_id = label_seq_id
                elif auth_id == end_res:
                    end_label_seq_id = label_seq_id
        
        logging.debug(f"Chain {chain_id} sequence mapping:")
        logging.debug(f"  Start residue {start_res} -> label_seq_id {start_label_seq_id}")
        logging.debug(f"  End residue {end_res} -> label_seq_id {end_label_seq_id}")
        
        if start_label_seq_id is not None and end_label_seq_id is not None:
            # Get sequence indices from label_seq_ids
            start_idx = None
            end_idx = None
            
            for (cid, label_id), seq_idx in self.label_seq_id_to_full_seq_idx.items():
                if cid == chain_id:
                    if label_id == start_label_seq_id:
                        start_idx = seq_idx
                    elif label_id == end_label_seq_id:
                        end_idx = seq_idx
            
            logging.debug(f"  Start label_seq_id {start_label_seq_id} -> index {start_idx}")
            logging.debug(f"  End label_seq_id {end_label_seq_id} -> index {end_idx}")
            logging.debug(f"  Full sequence length: {len(full_sequence)}")
            
            if start_idx is not None and end_idx is not None and start_idx < end_idx:
                # Extract the missing sequence
                missing_seq = full_sequence[start_idx + 1:end_idx]
                logging.debug(f"  Extracted sequence: {missing_seq}")
                return missing_seq
        
        # Fallback to X's if we can't get the actual sequence
        return 'X' * (end_res - start_res - 1)

    def load_structure(self) -> None:
        """Load and parse the structure file"""
        if self.pdb_file_path.endswith('.cif'):
            parser = MMCIFParser()
            self.cif_dict = MMCIF2Dict(self.pdb_file_path)
            self.structure = parser.get_structure('protein', self.pdb_file_path)
            self._parse_mmcif_sequences_and_mappings()
        else:
            parser = PDBParser()
            self.structure = parser.get_structure('protein', self.pdb_file_path)
            logging.warning("PDB files are less reliable for full sequence parsing. Missing loop sequences might be 'X's or require external fetching.")
            for model in self.structure:
                for chain in model:
                    seq_from_atom = ""
                    idx_counter = 0
                    for res in chain:
                        if is_aa(res):
                            try:
                                seq1_resname = seq1(res.resname)
                            except ValueError:
                                seq1_resname = 'X'
                            seq_from_atom += seq1_resname
                            self.observed_residues_by_label_seq_id[(chain.id, res.id[1])] = (res.id[1], res.id[2], res)
                            self.label_seq_id_to_full_seq_idx[(chain.id, res.id[1])] = idx_counter
                            self.original_flanking_residues[(chain.id, res.id[1], res.id[2])] = res
                            idx_counter += 1
                    self.full_chain_sequences[chain.id] = seq_from_atom

    def _parse_mmcif_sequences_and_mappings(self):
        """
        Extracts full chain sequences from mmCIF _entity_poly category
        and creates robust mappings for observed residues using Biopython's Structure object.
        """
        # Step 1: Extract full sequences from _entity_poly category
        if '_entity_poly.entity_id' in self.cif_dict:
            entity_ids = self.cif_dict['_entity_poly.entity_id']
            chain_ids_str_list = self.cif_dict['_entity_poly.pdbx_strand_id']
            one_letter_seqs = self.cif_dict['_entity_poly.pdbx_seq_one_letter_code']

            logging.debug("=== Entity Poly Data ===")
            for eid, cids_str, seq in zip(entity_ids, chain_ids_str_list, one_letter_seqs):
                logging.debug(f"Entity ID: {eid}")
                logging.debug(f"Chain IDs: {cids_str}")
                logging.debug(f"Sequence: {seq[:50]}...")  # Show first 50 chars

            # Populate self.full_chain_sequences
            for eid, cids_str, seq in zip(entity_ids, chain_ids_str_list, one_letter_seqs):
                full_seq_clean = seq.replace('\n', '').replace(' ', '')
                for cid in cids_str.replace(' ', '').split(','):
                    self.full_chain_sequences[cid] = full_seq_clean
                    logging.debug(f"Mapped chain {cid} to sequence: {full_seq_clean[:50]}...")
            logging.debug(f"Full chain sequences parsed from _entity_poly: {self.full_chain_sequences}")

        # Step 2: Create mapping from Biopython Structure
        logging.debug("\n=== Building Biopython Residue Map ===")
        
        # First, build the primary mapping from Biopython Structure
        temp_biopython_res_map = {}
        
        # Create a lookup dictionary for _atom_site records
        atom_site_lookup = {}
        if '_atom_site.label_seq_id' in self.cif_dict:
            logging.debug("\nProcessing _atom_site records:")
            for i in range(len(self.cif_dict['_atom_site.label_asym_id'])):
                auth_asym_id = self.cif_dict['_atom_site.auth_asym_id'][i]
                auth_seq_id = self.cif_dict['_atom_site.auth_seq_id'][i]
                auth_ins_code = self.cif_dict['_atom_site.pdbx_PDB_ins_code'][i]
                label_seq_id = self.cif_dict['_atom_site.label_seq_id'][i]
                
                # Skip if auth_seq_id is not a valid number
                try:
                    auth_seq_id_int = int(auth_seq_id)
                except ValueError:
                    logging.debug(f"Skipping invalid auth_seq_id: {auth_seq_id}")
                    continue
                
                # Skip if label_seq_id is not a valid number
                try:
                    label_seq_id_int = int(label_seq_id)
                except ValueError:
                    logging.debug(f"Skipping invalid label_seq_id: {label_seq_id}")
                    continue
                
                # Create a key that matches Biopython's residue ID format
                key = (auth_asym_id, auth_seq_id_int, auth_ins_code if auth_ins_code not in ['?', '.', ' '] else ' ')
                atom_site_lookup[key] = label_seq_id_int
                # Store the auth_id to label_seq_id mapping
                self.auth_id_to_label_seq_id[(auth_asym_id, auth_seq_id_int)] = label_seq_id_int
                logging.debug(f"Added to auth_id_to_label_seq_id: ({auth_asym_id}, {auth_seq_id_int}) -> {label_seq_id_int}")
        
        logging.debug("\nVerifying auth_id_to_label_seq_id mapping:")
        for (chain_id, auth_id), label_id in self.auth_id_to_label_seq_id.items():
            logging.debug(f"Chain {chain_id}, Auth ID {auth_id} -> Label ID {label_id}")
        
        for model in self.structure:
            for chain in model:
                chain_id = chain.id
                logging.debug(f"\nProcessing chain {chain_id}:")
                
                # Get all residues in the chain
                residues = list(chain.get_residues())
                residues.sort(key=lambda r: r.id[1])  # Sort by residue number
                
                for res in residues:
                    if not is_aa(res):
                        continue
                        
                    # Store with Biopython's normalized keys
                    key = (chain_id, res.id[1], res.id[2])
                    temp_biopython_res_map[key] = res
                    self.original_flanking_residues[key] = res
                    logging.debug(f"  Added residue: Chain {chain_id}, ResID {res.id}, ResName {res.resname}")
                    
                    # Get label_seq_id from our lookup dictionary
                    if key in atom_site_lookup:
                        label_seq_id = atom_site_lookup[key]
                        self.observed_residues_by_label_seq_id[(chain_id, label_seq_id)] = (res.id[1], res.id[2], res)
                        
                        # Map to full sequence index if possible
                        if chain_id in self.full_chain_sequences:
                            full_seq = self.full_chain_sequences[chain_id]
                            if label_seq_id - 1 >= 0 and label_seq_id - 1 < len(full_seq):
                                self.label_seq_id_to_full_seq_idx[(chain_id, label_seq_id)] = label_seq_id - 1
                                logging.debug(f"    Mapped to sequence index {label_seq_id - 1}")

        logging.debug("\n=== Summary ===")
        logging.debug(f"Total residues mapped: {len(temp_biopython_res_map)}")
        logging.debug(f"Observed residues by Label Seq ID: {len(self.observed_residues_by_label_seq_id)}")
        logging.debug(f"Label Seq ID to Full Seq Index mapping: {len(self.label_seq_id_to_full_seq_idx)}")
        logging.debug(f"Auth ID to Label Seq ID mapping: {len(self.auth_id_to_label_seq_id)}")
        
        # Debug output for sequence mapping
        for chain_id, seq in self.full_chain_sequences.items():
            logging.debug(f"\nChain {chain_id} sequence mapping:")
            logging.debug(f"Full sequence: {seq[:50]}...")
            mapped_indices = [idx for (cid, _), idx in self.label_seq_id_to_full_seq_idx.items() if cid == chain_id]
            if mapped_indices:
                logging.debug(f"Mapped indices: {sorted(mapped_indices)}")
            else:
                logging.debug("No indices mapped for this chain")

    def find_missing_regions(self) -> List[MissingRegion]:
        """Find all missing regions in the structure"""
        if not self.structure:
            raise ValueError("Structure not loaded. Call load_structure() first.")
        missing_regions = []
        
        for model in self.structure:
            for chain in model:
                chain_id = chain.id
                # Get all residues in the chain
                residues = list(chain.get_residues())
                if not residues:
                    continue
                
                # Sort residues by their ID number
                residues.sort(key=lambda r: r.id[1])
                
                # Find gaps between consecutive residues
                for i in range(len(residues) - 1):
                    curr_res = residues[i]
                    next_res = residues[i + 1]
                    
                    # Skip if not standard amino acids
                    if not is_aa(curr_res) or not is_aa(next_res):
                        continue
                    
                    # Check if there's a gap
                    curr_id = curr_res.id[1]
                    next_id = next_res.id[1]
                    if next_id - curr_id > 1:
                        # Found a gap
                        start_res = curr_id
                        end_res = next_id
                        length = end_res - start_res - 1
                        
                        # Get coordinates
                        try:
                            start_coords = np.array(curr_res['CA'].get_coord())
                        except KeyError:
                            start_coords = np.array(list(curr_res.get_atoms())[0].get_coord())
                        try:
                            end_coords = np.array(next_res['CA'].get_coord())
                        except KeyError:
                            end_coords = np.array(list(next_res.get_atoms())[0].get_coord())
                        
                        # Get sequence if available, otherwise use X's
                        full_chain_sequence = self.full_chain_sequences.get(chain_id, '')
                        missing_sequence = self._get_missing_sequence(chain_id, start_res, end_res)
                        
                        missing_regions.append(MissingRegion(
                            chain_id=chain_id,
                            start_res=start_res,
                            end_res=end_res,
                            length=length,
                            start_coords=start_coords,
                            end_coords=end_coords,
                            start_res_name=curr_res.resname,
                            end_res_name=next_res.resname,
                            missing_sequence=missing_sequence,
                            full_chain_sequence=full_chain_sequence,
                            start_icode=curr_res.id[2],
                            end_icode=next_res.id[2]
                        ))
                        logging.debug(f"Found missing region in chain {chain_id}: {start_res} to {end_res} (length: {length})")
                        logging.debug(f"Missing sequence: {missing_sequence}")
        
        self.missing_regions = missing_regions
        return missing_regions

    def _calculate_quality_score(self, coords: List[np.ndarray], region: MissingRegion) -> float:
        """Calculate a quality score for a conformation"""
        # Check distance from start/end points
        start_dist = np.linalg.norm(coords[0] - region.start_coords)
        end_dist = np.linalg.norm(coords[-1] - region.end_coords)
        
        # Check smoothness (angle between consecutive segments)
        angles = []
        for i in range(len(coords) - 2):
            v1 = coords[i+1] - coords[i]
            v2 = coords[i+2] - coords[i+1]
            angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            angles.append(angle)
        
        # Check bond lengths
        bond_lengths = [np.linalg.norm(coords[i+1] - coords[i]) for i in range(len(coords)-1)]
        
        # Calculate score (higher is better)
        score = 1.0
        score *= np.exp(-start_dist - end_dist)  # Penalize large distances from endpoints
        score *= np.exp(-np.std(angles))  # Penalize irregular angles
        score *= np.exp(-np.std(bond_lengths))  # Penalize irregular bond lengths
        
        return float(score)

    def generate_conformations(self, region: MissingRegion, num_conformations: int = 5) -> List[Dict]:
        """Generate multiple conformations for a missing region - to be implemented by subclasses"""
        raise NotImplementedError("Subclasses must implement generate_conformations")

    def save_results(self, output_dir: str) -> None:
        """Save predicted structures and metadata using the configured result writer"""
        if not hasattr(self, 'missing_regions'):
            raise ValueError("No missing regions found. Call find_missing_regions() first.")
            
        conformations = []
        for region in self.missing_regions:
            region_conformations = self.generate_conformations(region)
            conformations.extend(region_conformations)
            
        self.result_writer.write_results(output_dir, self.pdb_file_path, self.missing_regions, conformations)

    def _write_conformation(self, output_file: Path, region: MissingRegion, conf: Dict) -> None:
        """Write a conformation to a PDB file"""
        with open(output_file, 'w') as f:
            f.write(f"REMARK  Quality score: {conf['quality_score']:.3f}\n")
            for i, coord in enumerate(conf['coordinates']):
                res_num = region.start_res + i + 1
                f.write(f"ATOM  {i+1:5d}  CA  ALA {region.chain_id}{res_num:4d}    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}\n") 