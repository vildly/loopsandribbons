import sys
from pathlib import Path

# Add the parent directory to sys.path
sys.path.append(str(Path(__file__).parent.parent))

from loops_api.loop_predictor import LoopPredictor
import requests
import tempfile
import os
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser

sys.path.insert(0, "/opt/homebrew/Cellar/modeller/10.7/modlib")
import modeller
print(modeller.__version__)

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
        predictor = LoopPredictor(pdb_file)
        
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
                conformations = predictor.generate_conformations(region, num_conformations=2, use_modeller=True)
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
                
            except ImportError as e:
                print(f"Modeller not available: {e}")
                print("Skipping Modeller tests")
                break
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

def test_loop_prediction():
    # Get the structure file
    pdb_id = "3IDP"
    pdb_file = get_pdb_file(pdb_id)
    
    try:
        # Initialize predictor
        predictor = LoopPredictor(pdb_file)
        
        # Load structure and find missing regions
        predictor.load_structure()
        
        # Debug: Print mmCIF header data
        print("\nmmCIF Header Data:")
        if hasattr(predictor.structure, 'header'):
            header = predictor.structure.header
            if isinstance(header, dict):
                print("Available keys:", list(header.keys()))
                if '_entity_poly_seq' in header:
                    print("\n_entity_poly_seq records:")
                    for record in header['_entity_poly_seq']:
                        print(f"Chain {record.get('pdbx_strand_id')}: {record.get('mon_id', [])}")
        
        # Debug: Print full chain sequences
        print("\nFull chain sequences:")
        for chain_id, sequence in predictor.full_chain_sequences.items():
            print(f"Chain {chain_id}: {sequence}")
        
        missing_regions = predictor.find_missing_regions()
        
        print(f"\nFound {len(missing_regions)} missing regions in {pdb_id}:")
        for region in missing_regions:
            print(f"\nChain {region.chain_id}: {region.start_res} to {region.end_res} (length: {region.length})")
            print(f"  Start residue: {region.start_res_name}")
            print(f"  End residue: {region.end_res_name}")
            print(f"  Missing sequence: {region.missing_sequence}")
            print(f"  Full chain sequence length: {len(region.full_chain_sequence)}")
            
            # Verify sequence properties
            assert len(region.missing_sequence) == region.length, \
                f"Missing sequence length ({len(region.missing_sequence)}) doesn't match gap length ({region.length})"
            
            # Verify coordinates are present
            assert region.start_coords is not None and region.end_coords is not None, "Missing coordinates"
            assert region.start_coords.shape == (3,) and region.end_coords.shape == (3,), "Coordinates should be 3D vectors"
            
            # Verify the full chain sequence exists
            assert region.full_chain_sequence, "Full chain sequence should exist"
            
            # Verify the region is within the full chain sequence length
            assert region.start_res < region.end_res, "Start residue should be before end residue"
            assert region.length > 0, "Region length should be positive"
        
        # Generate and save predictions
        output_dir = Path("predictions")
        predictor.save_results(output_dir)
        
        print(f"\nPredictions saved to {output_dir}/")
        print("\nGenerated files:")
        for file in output_dir.glob("*.pdb"):
            print(f"  {file.name}")
        
        # Print metadata
        metadata_file = output_dir / "metadata.json"
        if metadata_file.exists():
            print("\nMetadata summary:")
            with open(metadata_file) as f:
                import json
                metadata = json.load(f)
                print(f"  PDB ID: {metadata['pdb_file_path']}")
                print(f"  Number of regions: {len(metadata['missing_regions'])}")
    
    except Exception as e:
        print(f"Error during testing: {e}")
        raise  # Re-raise the exception to see the full traceback

if __name__ == "__main__":
    test_loop_prediction()
    test_modeller_predictions()  # Add the new test 