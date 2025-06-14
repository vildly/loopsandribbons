import requests
from Bio import PDB
import tempfile
import os
from typing import List, Dict, Tuple
import subprocess
import json

def download_pdb(pdb_id: str) -> str:
    """Download PDB file and return the path to the downloaded file"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    
    # Save to temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')
    temp_file.write(response.content)
    temp_file.close()
    
    return temp_file.name

def find_missing_regions(pdb_file: str) -> List[Dict]:
    """Find missing regions in the structure"""
    parser = PDB.PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    
    missing_regions = []
    
    for model in structure:
        for chain in model:
            prev_res = None
            for residue in chain:
                if prev_res and residue.id[1] - prev_res.id[1] > 1:
                    # Found a gap
                    missing_regions.append({
                        'chain_id': chain.id,
                        'start_res': prev_res.id[1],
                        'end_res': residue.id[1],
                        'length': residue.id[1] - prev_res.id[1] - 1
                    })
                prev_res = residue
    
    return missing_regions

def prepare_rarefold_input(pdb_file: str, missing_region: Dict) -> str:
    """Prepare input for RareFold"""
    # Create a temporary directory for RareFold input
    temp_dir = tempfile.mkdtemp()
    
    # Copy the PDB file to the temp directory
    input_pdb = os.path.join(temp_dir, 'input.pdb')
    with open(pdb_file, 'r') as f:
        with open(input_pdb, 'w') as out:
            out.write(f.read())
    
    # Create a JSON file with the missing region information
    region_info = {
        'chain_id': missing_region['chain_id'],
        'start_res': missing_region['start_res'],
        'end_res': missing_region['end_res']
    }
    
    with open(os.path.join(temp_dir, 'region_info.json'), 'w') as f:
        json.dump(region_info, f)
    
    return temp_dir

def run_rarefold(input_dir: str) -> str:
    """Run RareFold prediction"""
    # This is a placeholder for the actual RareFold command
    # You'll need to modify this based on how RareFold is installed and used
    cmd = f"rarefold --input_dir {input_dir} --output_dir {input_dir}/output"
    subprocess.run(cmd, shell=True, check=True)
    
    return os.path.join(input_dir, 'output', 'predicted.pdb')

def main():
    # Download the PDB file
    pdb_id = "3IDP"
    pdb_file = download_pdb(pdb_id)
    
    try:
        # Find missing regions
        missing_regions = find_missing_regions(pdb_file)
        print(f"Found {len(missing_regions)} missing regions:")
        for region in missing_regions:
            print(f"Chain {region['chain_id']}: {region['start_res']} to {region['end_res']} (length: {region['length']})")
        
        # Process each missing region
        for region in missing_regions:
            print(f"\nProcessing region in chain {region['chain_id']}...")
            
            # Prepare input for RareFold
            input_dir = prepare_rarefold_input(pdb_file, region)
            
            # Run RareFold
            output_file = run_rarefold(input_dir)
            
            print(f"Prediction saved to: {output_file}")
            
    finally:
        # Clean up
        os.unlink(pdb_file)

if __name__ == "__main__":
    main() 