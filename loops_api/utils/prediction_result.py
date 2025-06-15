from abc import ABC, abstractmethod
from pathlib import Path
import json
from typing import List, Dict, Optional
from Bio.PDB import PDBIO
from .missing_region import MissingRegion
import datetime

class ResultWriter(ABC):
    """Abstract base class for writing prediction results"""
    
    @abstractmethod
    def write_results(self, output_dir: str, pdb_file_path: str, region: MissingRegion, conformations: List[Dict]) -> Path:
        """Write prediction results for a single region to the specified output location"""
        pass
    
    @abstractmethod
    def write_metadata(self, output_dir: str, pdb_file_path: str, missing_regions: List[MissingRegion]) -> None:
        """Write overall metadata for all regions"""
        pass

class FileResultWriter(ResultWriter):
    """Concrete implementation of ResultWriter that writes to files"""
    
    def __init__(self):
        self.base_dir = Path('predictions')
        self.base_dir.mkdir(exist_ok=True)
    
    def _create_prediction_dir(self) -> Path:
        """Create a new prediction directory with timestamp"""
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        prediction_dir = self.base_dir / f"loop_prediction_{timestamp}"
        prediction_dir.mkdir(exist_ok=True, parents=True)
        return prediction_dir
    
    def write_results(self, output_dir: str, pdb_file_path: str, region: MissingRegion, conformations: List[Dict]) -> Path:
        """Write predicted structures for a single region to files"""
        # Use the provided output directory and ensure it exists
        print(output_dir)
        prediction_summary_file = Path("prediction_summary.json")
     
        # Save prediction summary for this region
        summary = {
            'timestamp': datetime.datetime.now().strftime("%Y%m%d_%H%M%S"),
            'region_metadata': {
                'chain_id': region.chain_id,
                'start_res': region.start_res,
                'end_res': region.end_res,
                'length': region.length,
                'missing_sequence': region.missing_sequence,
                'full_chain_sequence': region.full_chain_sequence
            },
            'num_models': len(conformations),
            'models': [
                {
                    'model_number': conf['conformation_id'],
                    'quality_score': conf['quality_score'],
                    # Add any additional scores if present
                    **{k: v for k, v in conf.items() if k.endswith('_score') and k != 'quality_score'},
                    'model_file': conf.get('model_file', '')
                }
                for conf in conformations
            ]
        }
        # save summary to file        
        with open(prediction_summary_file, 'w+') as f:
            json.dump(summary, f, indent=2)

        return prediction_summary_file

    def write_metadata(self, output_dir: str, pdb_file_path: str, missing_regions: List[MissingRegion]) -> None:
        """Write overall metadata for all regions"""
        output_dir = Path(output_dir)
        
        # Write standard metadata
        metadata = {
            'pdb_file_path': pdb_file_path,
            'missing_regions': [
                {
                    'chain_id': region.chain_id,
                    'start_res': region.start_res,
                    'end_res': region.end_res,
                    'length': region.length,
                    'start_res_name': region.start_res_name,
                    'end_res_name': region.end_res_name,
                    'missing_sequence': region.missing_sequence,
                    'full_chain_sequence': region.full_chain_sequence
                }
                for region in missing_regions
            ]
        }
        
        # Save metadata
        with open(output_dir / 'metadata.json', 'w') as f:
            json.dump(metadata, f, indent=2)

    def _write_conformation(self, output_file: Path, region: MissingRegion, conf: Dict) -> None:
        """Write a conformation to a PDB file"""
        with open(output_file, 'w') as f:
            # Write all available scores as remarks
            for key, value in conf.items():
                if key.endswith('_score'):
                    f.write(f"REMARK  {key.replace('_', ' ').title()}: {value:.3f}\n")
            
            # Write coordinates
            for i, coord in enumerate(conf['coordinates']):
                res_num = region.start_res + i + 1
                f.write(f"ATOM  {i+1:5d}  CA  ALA {region.chain_id}{res_num:4d}    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}\n") 