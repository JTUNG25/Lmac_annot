#!/usr/bin/env python3
"""
Rename gene IDs in GFF3 file to follow standard nomenclature.
Usage: python3 rename_gff_ids.py input.gff3 output.gff3 --prefix Lmac_D5
"""

import argparse
import re
from collections import defaultdict

def rename_gff_ids(input_gff, output_gff, prefix, zero_padding=5):
    """
    Rename GFF3 IDs following pattern: prefix_XXXXX
    """
    
    # Map old IDs to new IDs
    id_mapping = {}
    gene_counter = 0
    
    with open(input_gff, 'r') as f_in, open(output_gff, 'w') as f_out:
        for line in f_in:
            # Keep header lines unchanged
            if line.startswith('#'):
                f_out.write(line)
                continue
            
            # Skip empty lines
            if not line.strip():
                f_out.write(line)
                continue
            
            fields = line.rstrip('\n').split('\t')
            
            if len(fields) < 9:
                f_out.write(line)
                continue
            
            attributes = fields[8]
            
            # Extract ID and Parent information
            id_match = re.search(r'ID=([^;]+)', attributes)
            parent_match = re.search(r'Parent=([^;]+)', attributes)
            
            if id_match:
                old_id = id_match.group(1)
                
                # If it's a gene feature, create new ID with counter
                if fields[2] == 'gene':
                    gene_counter += 1
                    new_id = f"{prefix}_{str(gene_counter).zfill(zero_padding)}"
                    id_mapping[old_id] = new_id
                
                # Use mapped ID if available, otherwise keep as is
                new_id = id_mapping.get(old_id, old_id)
                
                # Update ID in attributes
                attributes = re.sub(
                    rf'ID={re.escape(old_id)}',
                    f'ID={new_id}',
                    attributes
                )
                
                # Update Parent references if they exist in mapping
                if parent_match:
                    old_parent = parent_match.group(1)
                    if old_parent in id_mapping:
                        new_parent = id_mapping[old_parent]
                        attributes = re.sub(
                            rf'Parent={re.escape(old_parent)}',
                            f'Parent={new_parent}',
                            attributes
                        )
            
            fields[8] = attributes
            f_out.write('\t'.join(fields) + '\n')
    
    print(f"✓ Renamed {gene_counter} genes")
    print(f"✓ Output written to: {output_gff}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Rename gene IDs in GFF3 file'
    )
    parser.add_argument('input', help='Input GFF3 file')
    parser.add_argument('output', help='Output GFF3 file')
    parser.add_argument(
        '--prefix',
        default='Lmac_D5',
        help='Prefix for new gene IDs (default: Lmac_D5)'
    )
    parser.add_argument(
        '--padding',
        type=int,
        default=5,
        help='Zero-padding for gene numbers (default: 5, gives 00001)'
    )
    
    args = parser.parse_args()
    rename_gff_ids(args.input, args.output, args.prefix, args.padding)