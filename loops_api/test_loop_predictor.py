import sys
from pathlib import Path

# Add the parent directory to sys.path
sys.path.append(str(Path(__file__).parent.parent))

from utils.base_loop_predictor import LoopPredictor, MissingRegion
from utils.simple_predictor import SimplePredictor
from utils.modeller_predictor import ModellerPredictor, MODELLER_AVAILABLE

import requests
import tempfile
import os
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

# Add Modeller path if available
if MODELLER_AVAILABLE:
    sys.path.insert(0, "/opt/homebrew/Cellar/modeller/10.7/modlib")
    import modeller
    print(f"Using Modeller version {modeller.__version__}")

def get_pdb_file(pdb_id: str) -> str:
    """Get PDB file path, downloading if necessary"""
    # Check if we already have the mmCIF file
    cif_file = Path(f"{pdb_id}.cif")
    if cif_file.exists():
        return str(cif_file)
    
    # Try mmCIF first
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    response = requests.get(url)
    
    if response.status_code == 200:
        # Save to file
        with open(cif_file, 'wb') as f:
            f.write(response.content)
        return str(cif_file)
    
    # Fallback to PDB
    pdb_file = Path(f"{pdb_id}.pdb")
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    
    # Save to file
    with open(pdb_file, 'wb') as f:
        f.write(response.content)
    
    return str(pdb_file)

def extract_sequences_from_cif(cif_path):
    mmcif_dict = MMCIF2Dict(cif_path)
    # Get mapping from entity_id to chain_id(s)
    entity_ids = mmcif_dict['_entity_poly.entity_id']
    chain_ids = mmcif_dict['_entity_poly.pdbx_strand_id']
    seqs = mmcif_dict['_entity_poly.pdbx_seq_one_letter_code']
    entity_to_chains = {}
    for eid, cids in zip(entity_ids, chain_ids):
        for cid in cids.replace(' ', '').split(','):
            entity_to_chains[cid] = eid
    # Map chain_id to sequence
    chain_to_seq = {}
    for eid, cids, seq in zip(entity_ids, chain_ids, seqs):
        for cid in cids.replace(' ', '').split(','):
            chain_to_seq[cid] = seq.replace('\n', '').replace(' ', '')
    return chain_to_seq

def build_resnum_to_seqidx_map(chain, full_chain_sequence):
    # Build a list of (resnum, one_letter) from ATOM records
    pdb_residues = []
    for residue in chain:
        if is_aa(residue):
            resnum = residue.id[1]
            try:
                aa = seq1(residue.resname)
            except Exception:
                aa = 'X'
            pdb_residues.append((resnum, aa))
    # Now, align this to the full_chain_sequence
    seq = ''.join(aa for _, aa in pdb_residues)
    idx = full_chain_sequence.find(seq)
    if idx == -1:
        # fallback: try to align by first residue
        idx = 0
    mapping = {}
    for i, (resnum, _) in enumerate(pdb_residues):
        mapping[resnum] = idx + i
    return mapping

def build_full_resnum_to_seqidx_map(mmcif_dict, chain_id, entity_id):
    # Get all sequence numbers and mon_ids for this entity
    nums = [int(n) for eid, n in zip(mmcif_dict['_entity_poly_seq.entity_id'], mmcif_dict['_entity_poly_seq.num']) if eid == entity_id]
    mon_ids = [m for eid, m in zip(mmcif_dict['_entity_poly_seq.entity_id'], mmcif_dict['_entity_poly_seq.mon_id']) if eid == entity_id]
    # Map PDB residue numbers to sequence indices
    mapping = {}
    for seqidx, resnum in enumerate(nums):
        mapping[resnum] = seqidx
    return mapping

def test_modeller_predictions():
    """Test Modeller-based loop predictions and structure insertion"""
    if not MODELLER_AVAILABLE:
        print("Modeller is not available. Skipping Modeller tests.")
        return
        
    # Get the structure file
    pdb_id = "3IDP"  # Using 3IDP as it has a good example of a missing loop
    pdb_file = get_pdb_file(pdb_id)
    
    try:
        # Check for Modeller license key
        if not os.environ.get('MODLLER_KEY'):
            # Try to find the license key in the Modeller installation
            try:
                modeller_dir = os.path.dirname(modeller.__file__)
                license_file = os.path.join(modeller_dir, 'modlib', 'modeller.license')
                if os.path.exists(license_file):
                    with open(license_file, 'r') as f:
                        license_key = f.read().strip()
                        os.environ['MODLLER_KEY'] = license_key
                        print("Found Modeller license key in installation directory")
                else:
                    print("No Modeller license key found. Please set MODLLER_KEY environment variable.")
                    print("You can get a license key from https://salilab.org/modeller/registration.html")
                    return
            except ImportError:
                print("Modeller not installed. Please install it first.")
                return
        
        # Initialize predictor
        predictor = ModellerPredictor(pdb_file)
        
        # Load structure and find missing regions
        predictor.load_structure()
        missing_regions = predictor.find_missing_regions()
        
        print(f"\nTesting Modeller predictions for {pdb_id}")
        print(f"Found {len(missing_regions)} missing regions")
        
        # Test each missing region
        for i, region in enumerate(missing_regions):
            print(f"\nTesting region {i+1}:")
            print(f"Chain {region.chain_id}: {region.start_res} to {region.end_res} (length: {region.length})")
            print(f"Missing sequence: {region.missing_sequence}")
            
            # Generate predictions using Modeller
            try:
                conformations = predictor.generate_conformations(region, num_conformations=2)
                print(f"Generated {len(conformations)} Modeller conformations")
                
                # Verify each conformation
                for conf in conformations:
                    # Check that we have coordinates
                    assert 'coordinates' in conf, "Missing coordinates in conformation"
                    assert len(conf['coordinates']) == region.length, \
                        f"Number of coordinates ({len(conf['coordinates'])}) doesn't match region length ({region.length})"
                    
                    # Check that we have a complete structure
                    assert 'complete_structure' in conf, "Missing complete structure in conformation"
                    complete_structure = conf['complete_structure']
                    
                    # Verify the structure has the correct chain
                    assert region.chain_id in complete_structure[0], \
                        f"Chain {region.chain_id} not found in complete structure"
                    
                    # Get the chain and verify it has all residues
                    chain = complete_structure[0][region.chain_id]
                    residue_numbers = [res.id[1] for res in chain]
                    
                    # Check that the missing region is present
                    for res_num in range(region.start_res, region.end_res + 1):
                        assert res_num in residue_numbers, \
                            f"Residue {res_num} not found in complete structure"
                    
                    # Check that the loop coordinates match the sequence
                    loop_coords = conf['coordinates']
                    assert len(loop_coords) == len(region.missing_sequence), \
                        f"Number of loop coordinates ({len(loop_coords)}) doesn't match sequence length ({len(region.missing_sequence)})"
                    
                    print(f"Conformation {conf['conformation_id']} verified successfully")
                    print(f"Quality score: {conf['quality_score']:.3f}")
                
            except Exception as e:
                print(f"Error during Modeller prediction: {e}")
                raise
        
        # Save results
        output_dir = Path("predictions")
        predictor.save_results(output_dir)
        
        print(f"\nPredictions saved to {output_dir}/")
        print("\nGenerated files:")
        for file in output_dir.glob("*.pdb"):
            print(f"  {file.name}")
        
    except Exception as e:
        print(f"Error during testing: {e}")
        raise

def test_simple_predictions():
    """Test simple linear interpolation predictions"""
    # Get the structure file
    pdb_id = "3IDP"
    pdb_file = get_pdb_file(pdb_id)
    
    try:
        # Initialize predictor
        predictor = SimplePredictor(pdb_file)
        
        # Load structure and find missing regions
        predictor.load_structure()
        missing_regions = predictor.find_missing_regions()
        
        print(f"\nTesting simple predictions for {pdb_id}")
        print(f"Found {len(missing_regions)} missing regions")
        
        # Test each missing region
        for i, region in enumerate(missing_regions):
            print(f"\nTesting region {i+1}:")
            print(f"Chain {region.chain_id}: {region.start_res} to {region.end_res} (length: {region.length})")
            print(f"Missing sequence: {region.missing_sequence}")
            
            # Generate predictions
            conformations = predictor.generate_conformations(region, num_conformations=5)
            print(f"Generated {len(conformations)} simple conformations")
            
            # Verify each conformation
            for conf in conformations:
                # Check that we have coordinates
                assert 'coordinates' in conf, "Missing coordinates in conformation"
                assert len(conf['coordinates']) == region.length, \
                    f"Number of coordinates ({len(conf['coordinates'])}) doesn't match region length ({region.length})"
                
                # Check that we have a quality score
                assert 'quality_score' in conf, "Missing quality score in conformation"
                assert isinstance(conf['quality_score'], float), "Quality score should be a float"
                
                print(f"Conformation {conf['conformation_id']} verified successfully")
                print(f"Quality score: {conf['quality_score']:.3f}")
        
        # Save results
        output_dir = Path("predictions")
        predictor.save_results(output_dir)
        
        print(f"\nPredictions saved to {output_dir}/")
        print("\nGenerated files:")
        for file in output_dir.glob("*.pdb"):
            print(f"  {file.name}")
        
    except Exception as e:
        print(f"Error during testing: {e}")
        raise

if __name__ == "__main__":
    #test_simple_predictions()
    test_modeller_predictions() 