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

# Configure logging
def setup_logging(prediction_dir: Path) -> logging.Logger:
    """Set up logging to write to a file in the prediction directory"""
    log_file = prediction_dir / "prediction.log"
    logger = logging.getLogger("modeller_predictor")
    logger.setLevel(logging.DEBUG)
    
    # Create file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    
    # Create console handler with a higher log level
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.WARNING)
    
    # Create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return logger

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
        self.temp_dir = Path(tempfile.mkdtemp())
        self.logger = None  # Will be initialized in predict_loop

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
        
        # Create predictions directory if it doesn't exist
        predictions_dir = Path('predictions')
        predictions_dir.mkdir(exist_ok=True)
        
        # Create a subdirectory for this specific prediction
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        prediction_dir = predictions_dir / f"loop_prediction_{timestamp}"
        prediction_dir.mkdir(exist_ok=True, parents=True)  # Add parents=True to create all necessary parent directories
        
        # Set up logging for this prediction
        self.logger = setup_logging(prediction_dir)
        self.logger.info("=== Starting Modeller Loop Prediction ===")
        
        # Save current directory and change to prediction directory
        original_dir = os.getcwd()
        os.chdir(prediction_dir)
        
        try:
            # Update Modeller's environment to use our directory
            self.env.io.atom_files_directory = [str(prediction_dir), original_dir]
            
            # Create alignment file and get template ID
            alignment_file, template_pdb_id = self._create_alignment_file(region, region.full_chain_sequence)
            self.logger.info(f"Using template ID: {template_pdb_id}")
            
            # Get gap indices from the alignment file
            with open(alignment_file, 'r') as f:
                lines = f.readlines()
                template_seq = ''.join(line.strip() for line in lines if not line.startswith('>') and not ':' in line)
                gap_idxs = self._get_gap_indexes(template_seq)
                self.logger.info(f"Gap indices from alignment: {gap_idxs}")
                
                if not gap_idxs:
                    raise ValueError("Could not find gap indices in alignment file")
            
            self.logger.info("=== Creating Modeller Loop Model ===")
            
            # Create a custom model class that only refines the missing region
            class LoopOnlyModel(LoopModel):
                def select_loop_atoms(self):
                    return Selection(self.residue_range(f'{gap_idxs[0]}:A', 
                                                      f'{gap_idxs[1]}:A'))
            
            # Initialize the model
            self.logger.info("Initializing Modeller automodel...")
            a = LoopOnlyModel(self.env, 
                             alnfile=alignment_file,
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
            
            # Modeller output base name (will be tempdir/basename.D0000000X.pdb)
            a.basename = self.temp_dir / f"loop_model_{region.chain_id}_{region.start_res}_{region.end_res}"
            
            self.logger.debug(f"Running MODELLER for loop {region.chain_id}:{region.start_res}-{region.end_res}")
            
            # Create log file in prediction directory
            modeller_output_log = Path("modeller.log")  # Use relative path since we're already in prediction_dir
            self.logger.info(f"Creating Modeller output log at: {modeller_output_log}")
            with open(modeller_output_log, 'w+') as log_f:
                old_stdout = sys.stdout
                sys.stdout = log_f 
                try:
                    a.make() # Run Modeller
                except Exception as e:
                    self.logger.error(f"Modeller failed with error: {e}")
                    raise
                finally:
                    sys.stdout = old_stdout 
            self.logger.info(f"MODELLER output log saved to {modeller_output_log}")
            
            self.logger.info("Processing results...")
            
            # Get list of successfully built models
            ok_models = [x for x in a.outputs if x['failure'] is None]
            
            # Sort models by DOPE score
            ok_models.sort(key=lambda x: x['DOPE score'])
            
            results = []
            for i, model in enumerate(ok_models[:num_models]):
                model_pdb_path = Path(model['name'])
                model_number = i + 1
                
                self.logger.info(f"Processing MODELLER model {model_number} from file: {model_pdb_path}")
                
                # Read the coordinates of the missing region
                coords = self._extract_coordinates(str(model_pdb_path), region)
                
                # Parse the complete structure from Modeller's output
                model_parser = PDBParser()
                complete_structure = model_parser.get_structure(f'modeller_model_{model_number}', str(model_pdb_path))
                
                # Get quality scores directly from Modeller's output
                quality_score = model['DOPE score']
                ga341_score = model.get('GA341 score', 0.0)  # GA341 score might not be present

                results.append({
                    'coordinates': coords,
                    'complete_structure': complete_structure,
                    'quality_score': quality_score,
                    'ga341_score': ga341_score,
                    'conformation_id': model_number,
                    'model_file': str(model_pdb_path)
                })
                
            if not results:
                self.logger.warning("No MODELLER models were processed. This might mean Modeller failed, or file patterns are wrong.")

            # Save a summary file
            summary = {
                'timestamp': timestamp,
                'region_metadata': {
                    'chain_id': region.chain_id,
                    'start_res': region.start_res,
                    'end_res': region.end_res,
                    'length': region.length,
                    'missing_sequence': region.missing_sequence,
                    'full_chain_sequence': region.full_chain_sequence
                },
                'num_models': num_models,
                'models': [{
                    'model_number': r['conformation_id'],
                    'quality_score': r['quality_score'],
                    'ga341_score': r['ga341_score'],
                    'model_file': r['model_file']
                } for r in results]
            }
            
            with open("prediction_summary.json", 'w') as f:
                json.dump(summary, f, indent=2)
            
            return results
            
        finally:
            # Always return to original directory
            os.chdir(original_dir)
            # Close all handlers
            if self.logger:
                for handler in self.logger.handlers[:]:
                    handler.close()
                    self.logger.removeHandler(handler)
    
    def _extract_coordinates(self, model_file: str, region: MissingRegion) -> List[np.ndarray]:
        """Extract CA coordinates from the Modeller output file"""
        coords = []
        self.logger.debug(f"Extracting coordinates from {model_file} for region {region.chain_id}:{region.start_res}-{region.end_res}")
        
        try:
            with open(model_file) as f:
                for line in f:
                    if line.startswith('ATOM') and 'CA' in line:
                        parts = line.split()
                        if len(parts) >= 9:  # Ensure we have enough parts
                            try:
                                # Modeller uses sequential numbering starting from 1
                                # We need to map this to our target region
                                res_num = int(parts[5])
                                # Calculate the offset from the start of our region
                                offset = res_num - 1  # Modeller starts at 1
                                if 0 <= offset < region.length:  # Check if this residue is in our target region
                                    x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                                    coords.append(np.array([x, y, z]))
                                    self.logger.debug(f"Found CA atom for residue {res_num} (offset {offset}) at ({x}, {y}, {z})")
                            except (ValueError, IndexError) as e:
                                self.logger.warning(f"Error parsing line in PDB file: {line.strip()} - {e}")
                                continue
        except Exception as e:
            self.logger.error(f"Error reading PDB file {model_file}: {e}")
            raise
            
        if not coords:
            self.logger.warning(f"No coordinates found in {model_file} for region {region.chain_id}:{region.start_res}-{region.end_res}")
        else:
            self.logger.info(f"Extracted {len(coords)} coordinates from {model_file}")
            
        return coords
    
    def _parse_modeller_dope_score(self, log_file_path: Path, model_number: int) -> float:
        """Parses the DOPE score for a specific model from Modeller's log file."""
        dope_score = -99999.0  # Default low score if not found
        try:
            if not log_file_path.exists():
                self.logger.warning(f"Modeller log file not found at {log_file_path}")
                return dope_score

            with open(log_file_path, 'r') as f:
                log_content = f.read()
                
                # Modeller logs DOPE score in the V-file and sometimes in the main log for the final models.
                # A robust way is to check the V-file first.
                # Fix: Extract last 4 digits from model number for V-file name
                v_file_number = model_number % 10000  # Get last 4 digits
                score_v_file = log_file_path.parent / f"model.V9999{v_file_number:04d}"  # e.g., model.V99990002 for model 10002
                self.logger.debug(f"Looking for V-file at: {score_v_file}")
                if score_v_file.exists():
                    try:
                        with open(score_v_file, 'r') as sf:
                            for line in sf:
                                if "DOPE score:" in line:
                                    match = re.search(r"DOPE score:\s*([-+]?\d*\.\d+|\d+)", line)
                                    if match:
                                        dope_score = float(match.group(1))
                                        self.logger.debug(f"Parsed DOPE score {dope_score} for model {model_number} from V-file {score_v_file}")
                                        return dope_score
                    except Exception as e:
                        self.logger.warning(f"Error parsing DOPE score from V-file {score_v_file}: {e}")
                else:
                    self.logger.debug(f"V-file not found at {score_v_file}, falling back to main log")
                
                # Fallback to main log file if V-file not found or parsing failed
                # Search for specific model output lines in the log
                # Example log line: "filename: loop_model_B_597_614.D00000001.pdb ... DOPE score: -1234.567"
                
                # Regex for final assessment of a specific model number, covering D, B, BL, DL prefixes
                model_specific_pattern = r"filename:\s*" + re.escape(log_file_path.stem.split('.')[0]) + r"\.([DB]|BL|DL)\d{7}" + re.escape(str(model_number)) + r"\.pdb.*?DOPE score:\s*([-+]?\d*\.\d+|\d+)"
                
                match = re.search(model_specific_pattern, log_content, re.DOTALL)
                if match:
                    dope_score = float(match.group(2))  # Group 2 is the score
                    self.logger.debug(f"Parsed DOPE score {dope_score} for model {model_number} from main log")
                else:
                    # Fallback: Find the last overall DOPE score if model-specific not found
                    matches = list(re.finditer(r'>>DOPE score:\s+(-?\d+\.\d+)', log_content))
                    if matches:
                        dope_score = float(matches[-1].group(1))  # Get the last one
                        self.logger.debug(f"Parsed last overall DOPE score {dope_score} from log (not model specific)")
                    else:
                        self.logger.warning(f"No DOPE score found in log for model {model_number}")

        except Exception as e:
            self.logger.warning(f"Error reading Modeller log file {log_file_path} for model {model_number}: {e}")
        return dope_score 