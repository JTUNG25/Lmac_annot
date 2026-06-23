#!/usr/bin/env python3
"""
Generate functional annotation report from combined annotation data.

Input:
  - combined: TSV with merged eggNOG and DIAMOND annotations

Output:
  - report: Excel workbook with multiple annotation sheets
  - summary: Text file with annotation completeness statistics
"""

import csv
from collections import defaultdict, Counter
from pathlib import Path


def read_combined_annotations(combined_file):
    """Read combined annotation TSV."""
    annotations = []
    
    with open(combined_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            annotations.append(row)
    
    return annotations


def parse_comma_separated(value):
    """Parse comma-separated values (e.g., GO terms)."""
    if not value or value in ('', 'N/A'):
        return []
    return [v.strip() for v in value.split(',') if v.strip()]


def generate_statistics(annotations):
    """
    Generate comprehensive statistics about annotations.
    """
    total_proteins = len(annotations)
    
    with_eggnog = sum(1 for a in annotations if a.get('has_eggnog_annotation') == 'Yes')
    with_diamond = sum(1 for a in annotations if a.get('has_diamond_hit') == 'Yes')
    with_both = sum(1 for a in annotations 
                    if a.get('has_eggnog_annotation') == 'Yes' and a.get('has_diamond_hit') == 'Yes')
    
    # Count GO terms
    go_terms_all = []
    for a in annotations:
        go_terms_all.extend(parse_comma_separated(a.get('go_terms', '')))
    unique_go_terms = set(go_terms_all)
    proteins_with_go = sum(1 for a in annotations if a.get('go_terms', '').strip())
    
    # Count KEGG pathways
    kegg_pathways_all = []
    for a in annotations:
        kegg_pathways_all.extend(parse_comma_separated(a.get('kegg_pathway', '')))
    unique_kegg_pathways = set(kegg_pathways_all)
    proteins_with_kegg = sum(1 for a in annotations if a.get('kegg_pathway', '').strip())
    
    # Count COG categories
    cog_categories = Counter()
    for a in annotations:
        cog = a.get('cog_category', '').strip()
        if cog:
            for cat in cog:  # Each character is a category
                cog_categories[cat] += 1
    proteins_with_cog = sum(1 for a in annotations if a.get('cog_category', '').strip())
    
    # CAZy families
    cazy_all = []
    for a in annotations:
        cazy_all.extend(parse_comma_separated(a.get('cazy_family', '')))
    unique_cazy = set(cazy_all)
    proteins_with_cazy = sum(1 for a in annotations if a.get('cazy_family', '').strip())
    
    # EC numbers
    ec_all = []
    for a in annotations:
        ec_all.extend(parse_comma_separated(a.get('ec_number', '')))
    unique_ec = set(ec_all)
    proteins_with_ec = sum(1 for a in annotations if a.get('ec_number', '').strip())
    
    # Taxonomy level distribution
    tax_levels = Counter()
    for a in annotations:
        tax = a.get('tax_level', '').strip()
        if tax and tax != 'N/A':
            tax_levels[tax] += 1
    
    # Identity distribution from DIAMOND hits
    pident_values = []
    for a in annotations:
        pident = a.get('diamond_pident', '').strip()
        if pident:
            try:
                pident_values.append(float(pident))
            except ValueError:
                pass
    
    stats = {
        'total_proteins': total_proteins,
        'with_eggnog': with_eggnog,
        'with_diamond': with_diamond,
        'with_both': with_both,
        'with_go_terms': proteins_with_go,
        'unique_go_terms': len(unique_go_terms),
        'with_kegg_pathways': proteins_with_kegg,
        'unique_kegg_pathways': len(unique_kegg_pathways),
        'with_cog': proteins_with_cog,
        'cog_categories': dict(cog_categories),
        'with_cazy': proteins_with_cazy,
        'unique_cazy': len(unique_cazy),
        'with_ec': proteins_with_ec,
        'unique_ec_numbers': len(unique_ec),
        'tax_level_distribution': dict(tax_levels),
        'diamond_pident_values': pident_values,
        'avg_pident': sum(pident_values) / len(pident_values) if pident_values else 0,
    }
    
    return stats


def write_summary_report(stats, output_file):
    """Write text summary report."""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("TIBERIUS FUNCTIONAL ANNOTATION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("OVERALL STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total proteins annotated:               {stats['total_proteins']:>8}\n")
        f.write(f"  With eggNOG annotations:             {stats['with_eggnog']:>8} ({100*stats['with_eggnog']/stats['total_proteins']:.1f}%)\n")
        f.write(f"  With DIAMOND BLAST hits:             {stats['with_diamond']:>8} ({100*stats['with_diamond']/stats['total_proteins']:.1f}%)\n")
        f.write(f"  With both eggNOG and DIAMOND:        {stats['with_both']:>8} ({100*stats['with_both']/stats['total_proteins']:.1f}%)\n")
        
        f.write("\nFUNCTIONAL CATEGORIES\n")
        f.write("-" * 80 + "\n")
        f.write(f"GO Terms:              {stats['with_go_terms']:>6} proteins assigned, {stats['unique_go_terms']:>6} unique terms\n")
        f.write(f"KEGG Pathways:         {stats['with_kegg_pathways']:>6} proteins assigned, {stats['unique_kegg_pathways']:>6} unique pathways\n")
        f.write(f"COG Categories:        {stats['with_cog']:>6} proteins assigned, {len(stats['cog_categories']):>6} unique categories\n")
        f.write(f"CAZy Families:         {stats['with_cazy']:>6} proteins assigned, {stats['unique_cazy']:>6} unique families\n")
        f.write(f"EC Numbers:            {stats['with_ec']:>6} proteins assigned, {stats['unique_ec_numbers']:>6} unique numbers\n")
        
        if stats['cog_categories']:
            f.write("\n  COG Category Distribution:\n")
            for cat in sorted(stats['cog_categories'].keys()):
                count = stats['cog_categories'][cat]
                f.write(f"    {cat}: {count}\n")
        
        f.write("\nTAXONOMY LEVEL DISTRIBUTION (eggNOG)\n")
        f.write("-" * 80 + "\n")
        if stats['tax_level_distribution']:
            for level in sorted(stats['tax_level_distribution'].keys()):
                count = stats['tax_level_distribution'][level]
                f.write(f"  {level}: {count}\n")
        
        f.write("\nDIAMOND BLAST HOMOLOGY STATISTICS\n")
        f.write("-" * 80 + "\n")
        if stats['diamond_pident_values']:
            f.write(f"Average percent identity:   {stats['avg_pident']:.2f}%\n")
            f.write(f"Min percent identity:       {min(stats['diamond_pident_values']):.2f}%\n")
            f.write(f"Max percent identity:       {max(stats['diamond_pident_values']):.2f}%\n")
        else:
            f.write("No DIAMOND matches found.\n")
        
        f.write("\n" + "=" * 80 + "\n")


def create_excel_workbook(annotations, stats, output_file):
    """
    Create Excel workbook with multiple sheets.
    Requires openpyxl and pandas.
    """
    try:
        import pandas as pd
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment
        from openpyxl.utils.dataframe import dataframe_to_rows
    except ImportError:
        print("WARNING: pandas or openpyxl not available. Skipping Excel workbook creation.")
        print("Install with: pip install pandas openpyxl")
        return
    
    # Convert annotations to DataFrame
    df = pd.DataFrame(annotations)
    
    # Create Excel writer
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Sheet 1: Full annotation table
        df.to_excel(writer, sheet_name='Annotations', index=False)
        
        # Sheet 2: Summary statistics
        summary_data = {
            'Metric': [
                'Total Proteins',
                'With eggNOG',
                'With DIAMOND',
                'With both eggNOG and DIAMOND',
                'With GO terms',
                'With KEGG pathways',
                'With COG categories',
                'With CAZy families',
                'With EC numbers',
                'Unique GO terms',
                'Unique KEGG pathways',
                'Unique CAZy families',
                'Unique EC numbers',
                'Avg DIAMOND percent identity',
            ],
            'Count': [
                stats['total_proteins'],
                stats['with_eggnog'],
                stats['with_diamond'],
                stats['with_both'],
                stats['with_go_terms'],
                stats['with_kegg_pathways'],
                stats['with_cog'],
                stats['with_cazy'],
                stats['with_ec'],
                stats['unique_go_terms'],
                stats['unique_kegg_pathways'],
                stats['unique_cazy'],
                stats['unique_ec_numbers'],
                f"{stats['avg_pident']:.2f}%",
            ],
            'Percentage': [
                '100%',
                f"{100*stats['with_eggnog']/stats['total_proteins']:.1f}%",
                f"{100*stats['with_diamond']/stats['total_proteins']:.1f}%",
                f"{100*stats['with_both']/stats['total_proteins']:.1f}%",
                f"{100*stats['with_go_terms']/stats['total_proteins']:.1f}%",
                f"{100*stats['with_kegg_pathways']/stats['total_proteins']:.1f}%",
                f"{100*stats['with_cog']/stats['total_proteins']:.1f}%",
                f"{100*stats['with_cazy']/stats['total_proteins']:.1f}%",
                f"{100*stats['with_ec']/stats['total_proteins']:.1f}%",
                '-',
                '-',
                '-',
                '-',
                '-',
            ],
        }
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        # Sheet 3: COG distribution
        if stats['cog_categories']:
            cog_data = {
                'COG Category': list(stats['cog_categories'].keys()),
                'Count': list(stats['cog_categories'].values()),
            }
            cog_df = pd.DataFrame(cog_data)
            cog_df = cog_df.sort_values('Count', ascending=False)
            cog_df.to_excel(writer, sheet_name='COG_Distribution', index=False)
        
        # Sheet 4: Taxonomy level distribution
        if stats['tax_level_distribution']:
            tax_data = {
                'Taxonomy Level': list(stats['tax_level_distribution'].keys()),
                'Count': list(stats['tax_level_distribution'].values()),
            }
            tax_df = pd.DataFrame(tax_data)
            tax_df = tax_df.sort_values('Count', ascending=False)
            tax_df.to_excel(writer, sheet_name='Taxonomy_Distribution', index=False)
    
    print(f"Excel report written to {output_file}")


if __name__ == '__main__':
    import snakemake
    
    combined_file = snakemake.input.combined
    report_file = snakemake.output.report
    summary_file = snakemake.output.summary
    
    # Read annotations
    print(f"Reading annotations from {combined_file}...")
    annotations = read_combined_annotations(combined_file)
    print(f"  ✓ Loaded {len(annotations)} annotations")
    
    # Generate statistics
    print("Generating statistics...")
    stats = generate_statistics(annotations)
    
    # Write text summary
    print(f"Writing summary to {summary_file}...")
    write_summary_report(stats, summary_file)
    
    # Create Excel workbook
    print(f"Creating Excel report at {report_file}...")
    create_excel_workbook(annotations, stats, report_file)
    
    print("\n✓ Report generation complete!")