import numpy as np
import tempfile
import os
import re
from typing import List, Dict, Tuple
import sys
from pathlib import Path
import logging
from Bio.PDB import PDBParser, MMCIFParser, Structure, Model, Chain, Residue, Atom, PDBIO
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
import shutil
from modeller import *
from modeller.automodel import *
from modeller.automodel import LoopModel
from modeller.optimizers import molecular_dynamics
import datetime
import json

# Add the parent directory to sys.path
sys.path.append(str(Path(__file__).parent.parent.parent))

from .base_loop_predictor import LoopPredictor, MissingRegion

# Try to import Modeller, but don't fail if it's not available
MODELLER_AVAILABLE = False
try:
    import modeller
    from modeller import *
    from modeller.automodel import *
    from modeller.automodel import LoopModel
    from modeller.optimizers import molecular_dynamics
    
    # Check for Modeller license key
    if not os.environ.get('MODLLER_KEY'):
        # Try to find the license key in the Modeller installation
        modeller_dir = os.path.dirname(modeller.__file__)
        license_file = os.path.join(modeller_dir, 'modlib', 'modeller.license')
        if os.path.exists(license_file):
            with open(license_file, 'r') as f:
                license_key = f.read().strip()
                os.environ['MODLLER_KEY'] = license_key
                logging.info("Found Modeller license key in installation directory")
        else:
            logging.warning("No Modeller license key found. Please set MODLLER_KEY environment variable.")
            raise ImportError("Modeller license key not found")
    
    MODELLER_AVAILABLE = True
    logging.info(f"Modeller version {modeller.__version__} successfully imported")
except ImportError as e:
    logging.warning(f"Modeller import failed: {e}")
    logging.warning("Modeller-based predictions will not work.")

class ModellerPredictor(LoopPredictor):
    """Loop prediction using Modeller"""
    
    def __init__(self, pdb_file_path: str):
        super().__init__(pdb_file_path)
        if not MODELLER_AVAILABLE:
            raise ImportError("Modeller is not available. Please install it from https://salilab.org/modeller/")
        
        self.env = modeller.Environ()
        self.env.io.atom_files_directory = [os.path.dirname(os.path.abspath(pdb_file_path)), '.', '../atom_files']
        self.env.io.hetatm = True
        self.env.io.water = True
        self.env.libs.topology.read('${LIB}/top_heav.lib')
        self.env.libs.parameters.read('${LIB}/par.lib')
        modeller.log.verbose()
        self.temp_dir = tempfile.mkdtemp()

    def generate_conformations(self, region: MissingRegion, num_conformations: int = 2) -> List[Dict]:
        """Generate conformations using Modeller"""
        return self.predict_loop(region, num_conformations)

    def _create_alignment_file(self, region: MissingRegion, template_seq: str) -> Tuple[str, str]:
        """Create an alignment file for Modeller using the validated MissingRegion data"""
        alignment_file = os.path.join(self.temp_dir, 'alignment.ali')
        
        # Get template ID from filename
        template_pdb_id = os.path.splitext(os.path.basename(self.pdb_file_path))[0].upper()
        template_pdb_filename = f"{template_pdb_id}.pdb"
        template_pdb_path = os.path.join(self.temp_dir, template_pdb_filename)
        
        # Save the structure in PDB format
        io = PDBIO()
        io.set_structure(self.structure)
        io.save(template_pdb_path)
        print(f"Saved template structure as: {template_pdb_path}")
        
        print("\n=== DETAILED SEQUENCE DEBUGGING ===")
        print(f"Template PDB ID: {template_pdb_id}")
        print(f"Looking for residues in chain {region.chain_id}")
        print(f"Expected start residue: {region.start_res}")
        print(f"Expected end residue: {region.end_res}")
        print(f"Expected missing sequence: {region.missing_sequence}")
        print(f"Full chain sequence: {region.full_chain_sequence}")
        
        # Get the chain from the structure
        chain_obj_in_pdb = None
        for model in self.structure:
            for chain in model:
                if chain.id == region.chain_id:
                    chain_obj_in_pdb = chain
                    break
            if chain_obj_in_pdb:
                break
        
        if not chain_obj_in_pdb:
            raise ValueError(f"Chain {region.chain_id} not found in structure for MODELLER input.")
        
        # Get first and last residue numbers from PDB
        sorted_chain_residues = sorted([r for r in chain_obj_in_pdb.get_residues() if is_aa(r)],
                                     key=lambda r: (r.id[1], r.id[2]))
        
        if not sorted_chain_residues:
            raise ValueError(f"No amino acid residues found in chain {region.chain_id} for MODELLER input.")
            
        first_auth_seq_id_in_pdb = sorted_chain_residues[0].id[1]
        last_auth_seq_id_in_pdb = sorted_chain_residues[-1].id[1]
        
        print(f"\nFirst residue in PDB: {first_auth_seq_id_in_pdb}")
        print(f"Last residue in PDB: {last_auth_seq_id_in_pdb}")
        
        # Build template sequence from actual PDB structure
        template_sequence_for_ali = []
        current_pos = first_auth_seq_id_in_pdb
        
        for res in sorted_chain_residues:
            res_num = res.id[1]
            
            # Add gaps for any missing residues before this one
            while current_pos < res_num:
                template_sequence_for_ali.append('-')
                current_pos += 1
            
            # Add the actual residue
            if region.start_res < res_num < region.end_res:
                # This is the loop we're modeling - add gap
                template_sequence_for_ali.append('-')
            else:
                # This residue is observed in PDB
                try:
                    aa = seq1(res.resname)
                except:
                    aa = 'X'
                template_sequence_for_ali.append(aa)
            current_pos = res_num + 1
        
        # Add any remaining gaps up to the last residue
        while current_pos <= last_auth_seq_id_in_pdb:
            template_sequence_for_ali.append('-')
            current_pos += 1
        
        template_sequence_with_gaps = "".join(template_sequence_for_ali)
        
        print(f"\nTemplate sequence with gaps: {template_sequence_with_gaps}")
        print(f"Length: {len(template_sequence_with_gaps)}")
        
        # Write the alignment file
        with open(alignment_file, 'w') as f:
            # Template sequence - use actual PDB residue numbers
            # structureX:3IDP:449:B:720:B::::
            f.write(f'>P1;{template_pdb_id}\n')
            f.write(f'structureX:{template_pdb_id}:{first_auth_seq_id_in_pdb}:{region.chain_id}:{last_auth_seq_id_in_pdb}:{region.chain_id}::::\n')
            f.write(f"{template_sequence_with_gaps}*\n")
            
            # Target sequence - use the same label_seq_id range as the template
            # Get the label_seq_id range that corresponds to the PDB residues
            first_label_seq_id = self.auth_id_to_label_seq_id.get((region.chain_id, first_auth_seq_id_in_pdb))
            if first_label_seq_id is None:
                raise ValueError(f"Could not map auth_seq_id {first_auth_seq_id_in_pdb} to label_seq_id")
            
            # New code - use the same label_seq_id range for both template and model
            f.write(f'>P1;model\n')
            f.write(f'sequence:model:{first_label_seq_id}:{region.chain_id}:{first_label_seq_id + len(template_sequence_with_gaps) - 1}:{region.chain_id}::::\n')
            # Trim the model sequence to match the template range
            start_idx = self.label_seq_id_to_full_seq_idx.get((region.chain_id, first_label_seq_id))
            if start_idx is None:
                raise ValueError(f"Could not map label_seq_id {first_label_seq_id} to sequence index")
            model_sequence = region.full_chain_sequence[start_idx:start_idx + len(template_sequence_with_gaps)]
            f.write(f"{model_sequence}*\n")
        
        # Debug logging for final alignment
        print("\nFinal alignment file content:")
        with open(alignment_file, 'r') as f:
            print(f.read())
            
        # Additional debugging
        print("\nAlignment details:")
        print(f"Template start/end: {first_auth_seq_id_in_pdb}:{region.chain_id}:{last_auth_seq_id_in_pdb}:{region.chain_id}")
        print(f"Model start/end: 1:{region.chain_id}:{len(region.full_chain_sequence)}:{region.chain_id}")
        print(f"Template sequence length: {len(template_sequence_with_gaps)}")
        print(f"Model sequence length: {len(region.full_chain_sequence)}")
        
        return alignment_file, template_pdb_id  # Return both the alignment file path and template ID
        
    def _remove_tails(self, alignment_file: str) -> None:
        """Remove leading and trailing gaps from the alignment file"""
        pattern = r':\s*-?\d+(\.\d+)?'
        seq_region = False

        with open(alignment_file, 'r') as ali:
            lines = ali.readlines()

        seq_num = 0
        first_line_in_fasta_seq = True
        len_start_gap = None
        len_end_gap = None
        for i, line in enumerate(lines):
            if re.search(pattern, line):
                seq_region = True
                seq_num += 1
                continue

            elif seq_region and line.startswith('-') and seq_num == 1:
                new_line = line.strip('-')
                len_start_gap = len(line) - len(new_line)
                lines[i] = new_line

            elif seq_region and line.strip().endswith('-*') and seq_num == 1:
                new_line = line.strip().strip('-*') + '*'
                len_end_gap = len(line) - len(new_line) - 1
                lines[i] = new_line
                seq_region = False

            elif seq_region and seq_num == 2:
                if first_line_in_fasta_seq and len_start_gap:
                    lines[i] = line[len_start_gap:]
                    first_line_in_fasta_seq = False

                if line.strip().endswith('*') and len_end_gap:
                    lines[i] = line.strip()[:-len_end_gap - 1] + '*' + '\n'
                    seq_region = False

        with open(alignment_file, 'w') as file:
            file.writelines(lines)
            
    def _get_gap_indexes(self, seq: str) -> List[int]:
        """Get the start and end indices of the gap in the sequence"""
        print("\nFinding gap indices in sequence:")
        print(f"Sequence: {seq}")
        
        # Find the first gap
        start_idx = None
        end_idx = None
        
        for i, char in enumerate(seq):
            if char == '-' and start_idx is None:
                start_idx = i
                print(f"Found start gap at index {i}")
            elif char != '-' and start_idx is not None and end_idx is None:
                end_idx = i
                print(f"Found end gap at index {i}")
                break
        
        if start_idx is None or end_idx is None:
            print("No gaps found in sequence!")
            return []
            
        print(f"Gap indices: {[start_idx, end_idx]}")
        return [start_idx, end_idx]
        
    def predict_loop(self, region: MissingRegion, num_models: int = 5) -> List[Dict]:
        """Predict loop structure using Modeller"""
        if not MODELLER_AVAILABLE:
            raise ImportError("Modeller is not available. Please install it from https://salilab.org/modeller/")
        
        print("\n=== Starting Modeller Loop Prediction ===")
        
        # Create predictions directory if it doesn't exist
        predictions_dir = Path('predictions')
        predictions_dir.mkdir(exist_ok=True)
        
        # Create a subdirectory for this specific prediction
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        prediction_dir = predictions_dir / f"loop_prediction_{timestamp}"
        prediction_dir.mkdir(exist_ok=True)
        
        # Save current directory and change to prediction directory
        original_dir = os.getcwd()
        os.chdir(prediction_dir)
        
        try:
            # Update Modeller's environment to use our directory
            self.env.io.atom_files_directory = [str(prediction_dir), original_dir]
            
            # Create alignment file and get template ID
            alignment_file, template_pdb_id = self._create_alignment_file(region, region.full_chain_sequence)
            print(f"Using template ID: {template_pdb_id}")
            
            # Get gap indices from the alignment file
            with open(alignment_file, 'r') as f:
                lines = f.readlines()
                template_seq = ''.join(line.strip() for line in lines if not line.startswith('>') and not ':' in line)
                gap_idxs = self._get_gap_indexes(template_seq)
                print(f"Gap indices from alignment: {gap_idxs}")
                
                if not gap_idxs:
                    raise ValueError("Could not find gap indices in alignment file")
            
            print("\n=== Creating Modeller Loop Model ===")
            
            # Create a custom model class that only refines the missing region
            class LoopOnlyModel(LoopModel):
                def select_loop_atoms(self):
                    return Selection(self.residue_range(f'{gap_idxs[0]}:A', 
                                                      f'{gap_idxs[1]}:A'))
            
            # Initialize the model
            print("\nInitializing Modeller automodel...")
            a = LoopOnlyModel(self.env, 
                             alnfile=alignment_file,  # Now relative to prediction_dir
                             knowns=template_pdb_id,
                             sequence='model',
                             assess_methods=(assess.DOPE, assess.GA341))
            
            # Set up model parameters
            a.starting_model = 1
            a.ending_model = num_models
            
            # Set up loop refinement parameters
            a.loop.starting_model = 1
            a.loop.ending_model = num_models
            a.loop.md_level = refine.fast
            
            print("\nStarting model generation...")
            # Generate models
            a.make()
            
            print("\nProcessing results...")
            # Process the results
            results = []
            for i in range(1, num_models + 1):
                model_file = f"model.B{i:05d}.pdb"  # Modeller uses 5 digits
                if os.path.exists(model_file):
                    print(f"Processing model {i}...")
                    # Read the coordinates of the missing region
                    coords = self._extract_coordinates(model_file, region)
                    # Create a new structure with the predicted loop inserted
                    complete_structure = self._insert_loop_into_structure(model_file, region)
                    
                    # Save the complete structure
                    complete_file = f"complete_model_{i}.pdb"
                    io = PDBIO()
                    io.set_structure(complete_structure)
                    io.save(complete_file)
                    
                    results.append({
                        'coordinates': coords,
                        'quality_score': self._calculate_quality_score(coords, region),
                        'conformation_id': i,
                        'complete_structure': complete_structure,
                        'model_file': str(prediction_dir / model_file),
                        'complete_file': str(prediction_dir / complete_file)
                    })
                    print(f"Processed model {i}")
                else:
                    print(f"Warning: Model file {model_file} not found")
            
            # Save a summary file
            summary = {
                'timestamp': timestamp,
                'chain_id': region.chain_id,
                'start_res': region.start_res,
                'end_res': region.end_res,
                'missing_sequence': region.missing_sequence,
                'num_models': num_models,
                'models': [{
                    'model_number': r['conformation_id'],
                    'quality_score': r['quality_score'],
                    'model_file': os.path.basename(r['model_file']),
                    'complete_file': os.path.basename(r['complete_file'])
                } for r in results]
            }
            
            with open("prediction_summary.json", 'w') as f:
                json.dump(summary, f, indent=2)
            
            return results
            
        finally:
            # Always return to original directory
            os.chdir(original_dir)
    
    def _extract_coordinates(self, model_file: str, region: MissingRegion) -> List[np.ndarray]:
        """Extract CA coordinates from the Modeller output file"""
        coords = []
        with open(model_file) as f:
            for line in f:
                if line.startswith('ATOM') and 'CA' in line:
                    parts = line.split()
                    res_num = int(parts[5])
                    if region.start_res <= res_num <= region.end_res:
                        x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                        coords.append(np.array([x, y, z]))
        return coords
    
    def _insert_loop_into_structure(self, model_file: str, region: MissingRegion) -> Structure:
        """Insert the predicted loop into the original structure"""
        # Create a copy of the original structure
        new_structure = Structure.Structure('complete')
        new_model = Model.Model(0)
        new_structure.add(new_model)
        
        # Copy all chains from the original structure
        for chain in self.structure[0]:
            new_chain = Chain.Chain(chain.id)
            new_model.add(new_chain)
            
            # Copy all residues from the original chain
            for residue in chain:
                if not (chain.id == region.chain_id and 
                       region.start_res < residue.id[1] < region.end_res):
                    new_chain.add(residue.copy())
        
        # Add the predicted loop residues
        target_chain = new_model[region.chain_id]
        with open(model_file) as f:
            current_residue = None
            for line in f:
                if line.startswith('ATOM'):
                    parts = line.split()
                    res_num = int(parts[5])
                    if region.start_res <= res_num <= region.end_res:
                        if current_residue is None or current_residue.id[1] != res_num:
                            # Create new residue
                            current_residue = Residue.Residue(
                                (' ', res_num, ' '),
                                parts[3],  # resname
                                parts[4]   # segid
                            )
                            target_chain.add(current_residue)
                        
                        # Add atom to current residue
                        atom = Atom.Atom(
                            parts[2],  # name
                            np.array([float(parts[6]), float(parts[7]), float(parts[8])]),
                            float(parts[9]),  # occupancy
                            float(parts[10]),  # bfactor
                            parts[2],  # altloc
                            parts[2],  # fullname
                            int(parts[1])  # serial number
                        )
                        current_residue.add(atom)
        
        return new_structure
    
    def _calculate_quality_score(self, coords: List[np.ndarray], region: MissingRegion) -> float:
        """Calculate a quality score for the predicted loop"""
        if not coords:
            return 0.0
            
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