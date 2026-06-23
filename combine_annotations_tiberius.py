#!/usr/bin/env python3
"""
Combine eggNOG-mapper and DIAMOND BLAST results into a single functional annotation table.

Input:
  - eggnog: eggNOG-mapper .emapper.annotations file
  - diamond: DIAMOND BLAST TSV results
  - gff: GFF3 genome annotation file (for locus info)

Output:
  - combined: TSV with merged annotations and functional categories
"""

import csv
from pathlib import Path
from collections import defaultdict


def parse_eggnog_annotations(eggnog_file):
    """
    Parse eggNOG-mapper .emapper.annotations output.
    
    Expected columns (tab-separated):
    0: query_name
    1: seed_eggNOG_ortholog
    2: seed_ortholog_evalue
    3: seed_ortholog_score
    4: best_tax_level
    5: preferred_name
    6: GOs
    7: EC
    8: KEGG_ko
    9: KEGG_Pathway
    10: KEGG_Module
    11: KEGG_Reaction
    12: KEGG_rclass
    13: BRITE
    14: CAZy
    15: BiGG_Reaction
    16: COG_functional_category
    17: eggNOG_OGs
    18: bestOG|eggNOG|evalue|score
    19: COG_category (from consensus OGs)
    """
    eggnog_data = {}
    
    with open(eggnog_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        # Skip header comments
        for row in reader:
            if row[0].startswith('#'):
                continue
            if not row or len(row) < 3:
                continue
            
            query = row[0].strip()
            eggnog_data[query] = {
                'seed_ortholog': row[1].strip() if len(row) > 1 and row[1] else 'N/A',
                'seed_evalue': row[2].strip() if len(row) > 2 and row[2] else 'N/A',
                'seed_score': row[3].strip() if len(row) > 3 and row[3] else 'N/A',
                'tax_level': row[4].strip() if len(row) > 4 and row[4] else 'N/A',
                'preferred_name': row[5].strip() if len(row) > 5 and row[5] else 'N/A',
                'go_terms': row[6].strip() if len(row) > 6 and row[6] else '',
                'ec': row[7].strip() if len(row) > 7 and row[7] else '',
                'kegg_ko': row[8].strip() if len(row) > 8 and row[8] else '',
                'kegg_pathway': row[9].strip() if len(row) > 9 and row[9] else '',
                'kegg_module': row[10].strip() if len(row) > 10 and row[10] else '',
                'cazy': row[14].strip() if len(row) > 14 and row[14] else '',
                'cog_category': row[16].strip() if len(row) > 16 and row[16] else '',
            }
    
    return eggnog_data


def parse_diamond_results(diamond_file):
    """
    Parse DIAMOND BLAST results (TSV format).
    
    Expected columns (6 format with stitle):
    0: qseqid (query protein ID)
    1: sseqid (subject ID)
    2: pident (percent identity)
    3: length (alignment length)
    4: mismatch
    5: gapopen
    6-9: qstart, qend, sstart, send
    10: evalue
    11: bitscore
    12: stitle (subject description)
    """
    diamond_data = defaultdict(list)
    
    with open(diamond_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or len(row) < 13:
                continue
            
            query = row[0].strip()
            subject_desc = row[12].strip() if len(row) > 12 else 'Unknown'
            
            # Extract just the protein name from description (before first space or |)
            subject_name = subject_desc.split('|')[0] if '|' in subject_desc else subject_desc.split()[0]
            
            hit = {
                'subject_id': row[1].strip(),
                'subject_name': subject_name,
                'subject_desc': subject_desc,
                'pident': row[2].strip(),
                'length': row[3].strip(),
                'evalue': float(row[10].strip()),
                'bitscore': row[11].strip(),
            }
            
            diamond_data[query].append(hit)
    
    return diamond_data


def parse_gff(gff_file):
    """
    Extract locus information from GFF3 file.
    """
    locus_info = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'gene':
                continue
            
            attributes = fields[8]
            attrs = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    k, v = attr.split('=', 1)
                    attrs[k] = v
            
            gene_id = attrs.get('ID', '')
            if gene_id:
                locus_info[gene_id] = {
                    'seqname': fields[0],
                    'start': fields[3],
                    'end': fields[4],
                    'strand': fields[6],
                }
    
    return locus_info


def combine_and_output(eggnog_data, diamond_data, locus_info, output_file):
    """
    Combine eggNOG and DIAMOND data into a single TSV.
    """
    
    # Collect all query protein IDs
    all_queries = set(eggnog_data.keys()) | set(diamond_data.keys())
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(
            f,
            delimiter='\t',
            fieldnames=[
                'protein_id',
                'locus_tag',
                'seqname',
                'start',
                'end',
                'strand',
                'eggnog_ortholog',
                'eggnog_evalue',
                'eggnog_score',
                'tax_level',
                'preferred_name',
                'go_terms',
                'ec_number',
                'kegg_ko',
                'kegg_pathway',
                'kegg_module',
                'cazy_family',
                'cog_category',
                'best_diamond_hit',
                'diamond_subject_id',
                'diamond_pident',
                'diamond_evalue',
                'diamond_bitscore',
                'diamond_alignment_length',
                'diamond_description',
                'has_eggnog_annotation',
                'has_diamond_hit',
            ],
        )
        
        writer.writeheader()
        
        for query in sorted(all_queries):
            egg = eggnog_data.get(query, {})
            diam = diamond_data.get(query, [{}])[0] if diamond_data.get(query) else {}
            
            # Try to extract locus tag from protein ID (format: usually locus_tag_v1, etc.)
            locus_tag = query.rsplit('_', 1)[0] if '_' in query else query
            loc = locus_info.get(locus_tag, {})
            
            row = {
                'protein_id': query,
                'locus_tag': locus_tag,
                'seqname': loc.get('seqname', ''),
                'start': loc.get('start', ''),
                'end': loc.get('end', ''),
                'strand': loc.get('strand', ''),
                'eggnog_ortholog': egg.get('seed_ortholog', ''),
                'eggnog_evalue': egg.get('seed_evalue', ''),
                'eggnog_score': egg.get('seed_score', ''),
                'tax_level': egg.get('tax_level', ''),
                'preferred_name': egg.get('preferred_name', ''),
                'go_terms': egg.get('go_terms', ''),
                'ec_number': egg.get('ec', ''),
                'kegg_ko': egg.get('kegg_ko', ''),
                'kegg_pathway': egg.get('kegg_pathway', ''),
                'kegg_module': egg.get('kegg_module', ''),
                'cazy_family': egg.get('cazy', ''),
                'cog_category': egg.get('cog_category', ''),
                'best_diamond_hit': diam.get('subject_name', ''),
                'diamond_subject_id': diam.get('subject_id', ''),
                'diamond_pident': diam.get('pident', ''),
                'diamond_evalue': diam.get('evalue', ''),
                'diamond_bitscore': diam.get('bitscore', ''),
                'diamond_alignment_length': diam.get('length', ''),
                'diamond_description': diam.get('subject_desc', ''),
                'has_eggnog_annotation': 'Yes' if egg.get('seed_ortholog', '') != 'N/A' else 'No',
                'has_diamond_hit': 'Yes' if diam else 'No',
            }
            
            writer.writerow(row)


if __name__ == '__main__':
    import snakemake
    
    eggnog_file = snakemake.input.eggnog
    diamond_file = snakemake.input.diamond
    gff_file = snakemake.input.gff
    output_file = snakemake.output.combined
    log_file = snakemake.log[0]
    
    try:
        print(f"Parsing eggNOG results from {eggnog_file}...", file=open(log_file, 'a'))
        eggnog_data = parse_eggnog_annotations(eggnog_file)
        print(f"  ✓ Parsed {len(eggnog_data)} proteins", file=open(log_file, 'a'))
        
        print(f"\nParsing DIAMOND BLAST results from {diamond_file}...", file=open(log_file, 'a'))
        diamond_data = parse_diamond_results(diamond_file)
        print(f"  ✓ Parsed {len(diamond_data)} proteins with hits", file=open(log_file, 'a'))
        
        print(f"\nExtracting locus information from {gff_file}...", file=open(log_file, 'a'))
        locus_info = parse_gff(gff_file)
        print(f"  ✓ Extracted {len(locus_info)} loci", file=open(log_file, 'a'))
        
        print(f"\nCombining annotations and writing to {output_file}...", file=open(log_file, 'a'))
        combine_and_output(eggnog_data, diamond_data, locus_info, output_file)
        print(f"  ✓ Complete!", file=open(log_file, 'a'))
        
    except Exception as e:
        with open(log_file, 'a') as f:
            f.write(f"ERROR: {e}\n")
        raise