#!/usr/bin/env python3

import os
from pathlib import Path

# Config
GENOME = "Lmac_D5"
ASSEMBLY = "data/genomes/Lmac_D5.fasta"
GFF3 = "results/tiberius_evidence/tiberius_train.gff3"
OMARK_RESULT = "results/omark/Lmac_D5_evidence/Lmac_D5_evidence.omamer"
OMARK_DB = "data/omark/LUCA.h5"

OUTDIR = "results/contamination_screening"

rule target:
    input:
        f"{OUTDIR}/contamination_report.txt",
        f"{OUTDIR}/scaffold_summary.tsv",
        f"{OUTDIR}/species_breakdown.txt",

# ============================================================================
# STAGE 1: Parse OMAark output to identify contaminated proteins
# ============================================================================

rule parse_omark:
    """
    Parse OMAark .omamer file to extract species assignments.
    Identify proteins NOT from Leptosphaeria.
    """
    input:
        omark=OMARK_RESULT,
    output:
        contam_proteins=f"{OUTDIR}/contaminated_proteins_list.tsv",
        species_map=f"{OUTDIR}/protein_species_map.tsv",
    log:
        f"logs/parse_omark.log",
    script:
        """
        import pandas as pd
        from collections import defaultdict
        
        contaminated = {}  # protein_id → species
        all_proteins = {}   # protein_id → species
        target_species = {'Leptosphaeria maculans', 'Leptosphaeria', 'Lmac'}
        
        print("Parsing OMAark output...")
        
        # Parse OMAark output
        # Expecting tab-delimited with columns: query_id, hit_id, species, evalue, bitscore, ...
        # Adjust column indices based on your actual format
        
        with open(input.omark) as f:
            header = f.readline().strip()
            print(f"Header: {header}")
            
            for i, line in enumerate(f):
                fields = line.strip().split('\\t')
                if len(fields) < 3:
                    continue
                
                query_id = fields[0]
                # Try to find species info in available columns
                # Common formats: species is in column 2, or in the hit description
                species_info = fields[2] if len(fields) > 2 else "UNKNOWN"
                
                all_proteins[query_id] = species_info
                
                # Check if this is a non-target species
                is_target = any(target in species_info for target in target_species)
                if not is_target:
                    contaminated[query_id] = species_info
                
                if i < 5:  # Print first few lines for debugging
                    print(f"  Row {i}: {query_id} -> {species_info} (target={is_target})")
        
        # Write all protein-species map
        with open(output.species_map, 'w') as f:
            f.write("protein_id\\tspecies\\tcontamination_status\\n")
            for protein_id in sorted(all_proteins.keys()):
                species = all_proteins[protein_id]
                status = "CONTAMINATED" if protein_id in contaminated else "TARGET"
                f.write(f"{protein_id}\\t{species}\\t{status}\\n")
        
        # Write contaminated protein list
        with open(output.contam_proteins, 'w') as f:
            f.write("protein_id\\tspecies\\n")
            for protein_id in sorted(contaminated.keys()):
                f.write(f"{protein_id}\\t{contaminated[protein_id]}\\n")
        
        print(f"\\nTotal proteins: {len(all_proteins)}")
        print(f"Contaminated proteins: {len(contaminated)} ({100*len(contaminated)/len(all_proteins):.1f}%)")
        """

# ============================================================================
# STAGE 2: Load assembly and map proteins to scaffolds
# ============================================================================

rule map_proteins_to_scaffolds:
    """
    Parse GFF3 + FASTA to:
    1. Get scaffold sizes
    2. Map genes/proteins to scaffolds
    3. Count contaminated proteins per scaffold
    4. Generate detailed scaffold report
    """
    input:
        gff3=GFF3,
        assembly=ASSEMBLY,
        contam_proteins=f"{OUTDIR}/contaminated_proteins_list.tsv",
    output:
        scaffold_report=f"{OUTDIR}/contamination_report.txt",
        summary_tsv=f"{OUTDIR}/scaffold_summary.tsv",
    log:
        f"logs/map_proteins_to_scaffolds.log",
    script:
        """
        from Bio import SeqIO
        from collections import defaultdict
        import pandas as pd
        
        # Load scaffold sizes
        print("Loading scaffold sizes from FASTA...")
        scaffold_sizes = {}
        for record in SeqIO.parse(input.assembly, 'fasta'):
            scaffold_sizes[record.id] = len(record.seq)
        
        # Parse contaminated proteins
        print("Loading contaminated proteins...")
        contaminated = set()
        with open(input.contam_proteins) as f:
            next(f)  # Skip header
            for line in f:
                protein_id = line.strip().split('\\t')[0]
                contaminated.add(protein_id)
        
        print(f"Total contaminated proteins: {len(contaminated)}")
        
        # Parse GFF3 and map genes to scaffolds
        print("Parsing GFF3 and mapping proteins to scaffolds...")
        scaffold_proteins = defaultdict(list)
        gene_to_scaffold = {}
        
        with open(input.gff3) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\\t')
                if len(fields) < 9 or fields[2] not in ['gene', 'mRNA']:
                    continue
                
                scaffold = fields[0]
                attrs = fields[8]
                
                # Extract gene ID from attributes (GFF3 format: ID=...)
                gene_id = None
                for attr in attrs.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('ID=')[1].split('.')[0]  # Remove version
                        break
                
                if gene_id:
                    is_contam = gene_id in contaminated
                    scaffold_proteins[scaffold].append((gene_id, is_contam))
                    gene_to_scaffold[gene_id] = scaffold
        
        # Generate detailed report
        print("Generating contamination report...")
        report_lines = []
        summary_data = []
        
        report_lines.append("=" * 100)
        report_lines.append("CONTAMINATION ANALYSIS REPORT")
        report_lines.append("=" * 100)
        report_lines.append("")
        report_lines.append(f"Total scaffolds analyzed: {len(scaffold_proteins)}")
        report_lines.append(f"Total contaminated proteins: {len(contaminated)}")
        report_lines.append("")
        report_lines.append("=" * 100)
        report_lines.append("SCAFFOLD-BY-SCAFFOLD BREAKDOWN")
        report_lines.append("=" * 100)
        report_lines.append("")
        
        # Sort scaffolds by contamination percentage (descending)
        scaffold_stats = []
        for scaffold in sorted(scaffold_proteins.keys()):
            proteins = scaffold_proteins[scaffold]
            total = len(proteins)
            num_contam = sum(1 for _, is_contam in proteins if is_contam)
            pct_contam = (num_contam / total * 100) if total > 0 else 0
            size = scaffold_sizes.get(scaffold, 0)
            
            scaffold_stats.append({
                'scaffold': scaffold,
                'size_bp': size,
                'total_proteins': total,
                'contaminated_proteins': num_contam,
                'pct_contaminated': pct_contam,
            })
        
        # Sort by contamination percentage descending
        scaffold_stats_sorted = sorted(scaffold_stats, key=lambda x: -x['pct_contaminated'])
        
        # Print top contaminated scaffolds
        for stat in scaffold_stats_sorted:
            if stat['pct_contaminated'] > 0:  # Only print scaffolds with contamination
                report_lines.append(
                    f"Scaffold: {stat['scaffold']:<40} "
                    f"Size: {stat['size_bp']:>12,} bp   "
                    f"Proteins: {stat['total_proteins']:>4}   "
                    f"Contaminated: {stat['contaminated_proteins']:>4} ({stat['pct_contaminated']:>6.1f}%)"
                )
        
        report_lines.append("")
        report_lines.append("=" * 100)
        report_lines.append("CONTAMINATION CATEGORIES")
        report_lines.append("=" * 100)
        report_lines.append("")
        
        # Categorize scaffolds
        high_contam = [s for s in scaffold_stats_sorted if s['pct_contaminated'] >= 80]
        medium_contam = [s for s in scaffold_stats_sorted if 40 <= s['pct_contaminated'] < 80]
        low_contam = [s for s in scaffold_stats_sorted if 0 < s['pct_contaminated'] < 40]
        clean = [s for s in scaffold_stats_sorted if s['pct_contaminated'] == 0]
        
        report_lines.append(f"HIGH contamination (≥80%):    {len(high_contam)} scaffolds")
        for s in high_contam[:10]:  # Show top 10
            report_lines.append(
                f"  - {s['scaffold']:<40} {s['pct_contaminated']:>6.1f}% "
                f"({s['contaminated_proteins']}/{s['total_proteins']} proteins, {s['size_bp']:,} bp)"
            )
        if len(high_contam) > 10:
            report_lines.append(f"  ... and {len(high_contam) - 10} more")
        
        report_lines.append("")
        report_lines.append(f"MEDIUM contamination (40–80%): {len(medium_contam)} scaffolds")
        for s in medium_contam[:5]:  # Show top 5
            report_lines.append(
                f"  - {s['scaffold']:<40} {s['pct_contaminated']:>6.1f}% "
                f"({s['contaminated_proteins']}/{s['total_proteins']} proteins, {s['size_bp']:,} bp)"
            )
        if len(medium_contam) > 5:
            report_lines.append(f"  ... and {len(medium_contam) - 5} more")
        
        report_lines.append("")
        report_lines.append(f"LOW contamination (0–40%):    {len(low_contam)} scaffolds")
        if len(low_contam) > 0:
            report_lines.append(f"  (Likely annotation noise or assembly artifacts; generally safe to keep)")
        
        report_lines.append("")
        report_lines.append(f"CLEAN (0% contamination):     {len(clean)} scaffolds")
        
        report_lines.append("")
        report_lines.append("=" * 100)
        report_lines.append("RECOMMENDATIONS FOR MANUAL REVIEW")
        report_lines.append("=" * 100)
        report_lines.append("")
        report_lines.append("1. HIGH contamination (≥80%): Likely bacterial co-assembly. Strong candidates for removal.")
        report_lines.append("2. MEDIUM contamination (40–80%): Could be true contamination or assembly errors. Review before removing.")
        report_lines.append("3. LOW contamination (0–40%): Likely annotation noise. Generally safe to keep unless very small.")
        report_lines.append("4. Consider scaffold SIZE: Removing large contaminated scaffolds is higher priority than small ones.")
        report_lines.append("")
        report_lines.append("Next step: Review '02_scaffold_summary.tsv' for detailed breakdown.")
        report_lines.append("Then create a 'removal_list.txt' with scaffolds you want to remove.")
        
        # Write report
        with open(output.scaffold_report, 'w') as f:
            f.write('\\n'.join(report_lines))
        
        # Write TSV summary (for easy sorting/filtering in Excel or R)
        df = pd.DataFrame(scaffold_stats_sorted)
        df = df.sort_values('pct_contaminated', ascending=False)
        df.to_csv(output.summary_tsv, sep='\\t', index=False)
        
        print(f"✓ Report written to: {output.scaffold_report}")
        print(f"✓ Summary TSV written to: {output.summary_tsv}")
        print(f"\\nSummary:")
        print(f"  High contamination (≥80%):  {len(high_contam)} scaffolds")
        print(f"  Medium contamination:       {len(medium_contam)} scaffolds")
        print(f"  Low contamination:          {len(low_contam)} scaffolds")
        print(f"  Clean:                      {len(clean)} scaffolds")
        """

# ============================================================================
# STAGE 3: Species breakdown summary
# ============================================================================

rule species_breakdown:
    """
    Summarize which species are contaminating your assembly.
    """
    input:
        species_map=f"{OUTDIR}/protein_species_map.tsv",
    output:
        species_summary=f"{OUTDIR}/species_breakdown.txt",
    log:
        f"logs/species_breakdown.log",
    script:
        """
        from collections import defaultdict
        import pandas as pd
        
        species_count = defaultdict(int)
        status_count = defaultdict(int)
        
        with open(input.species_map) as f:
            next(f)  # Skip header
            for line in f:
                protein_id, species, status = line.strip().split('\\t')
                species_count[species] += 1
                status_count[status] += 1
        
        # Write summary
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("SPECIES COMPOSITION SUMMARY")
        report_lines.append("=" * 80)
        report_lines.append("")
        
        report_lines.append("Overall breakdown:")
        total = sum(status_count.values())
        for status, count in sorted(status_count.items(), key=lambda x: -x[1]):
            pct = 100 * count / total
            report_lines.append(f"  {status:<20} {count:>6} proteins ({pct:>5.1f}%)")
        
        report_lines.append("")
        report_lines.append("Contaminating species (top 20):")
        for species, count in sorted(species_count.items(), key=lambda x: -x[1])[:20]:
            pct = 100 * count / total
            is_target = 'TARGET' if 'Leptosphaeria' in species else 'CONTAMINANT'
            report_lines.append(f"  {species:<50} {count:>6} proteins ({pct:>5.1f}%) [{is_target}]")
        
        report_lines.append("")
        report_lines.append("=" * 80)
        
        with open(output.species_summary, 'w') as f:
            f.write('\\n'.join(report_lines))
        
        print('\\n'.join(report_lines))
        """
