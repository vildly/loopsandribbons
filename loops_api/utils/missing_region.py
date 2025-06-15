from dataclasses import dataclass
import numpy as np

@dataclass
class MissingRegion:
    chain_id: str
    start_res: int
    end_res: int
    length: int
    start_coords: np.ndarray
    end_coords: np.ndarray
    start_res_name: str
    end_res_name: str
    missing_sequence: str     # The actual 1-letter sequence of the missing residues
    full_chain_sequence: str  # The complete 1-letter sequence of the chain
    start_icode: str
    end_icode: str 